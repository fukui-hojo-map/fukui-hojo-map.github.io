"""Earth Engine を使わない雲の判定（Cloud Score+ の判断を学習した判定器）。s2open.py から使う。

正解: 公開済みの圃場内マップ px_<年>.png（Earth Engine の Cloud Score+ で晴れの田の画素に値、雲は 0）。
特徴: その画素の帯の値（AWS の sentinel-cogs と同じ L2A。反射率＝値/10000）だけから作る。前後の日・広い範囲は使わない
     （毎日の更新で、その日の画像だけで判定できるように。学習の試しで、前後の日の基準・雲からの距離・1km の雲の割合を抜いても精度は変わらなかった）。
判定器: scikit-learn の HistGradientBoostingClassifier。models/s2cloud_<版>.pkl と、説明の .json（特徴・しきい値・学習の数字・scikit-learn の版）。
田の画素で学習したので、田（筆ポリゴンの中）の画素にだけ使う（森・町・水の上は学習していない）。

  学習: python pipeline/s2cloud_train.py --data <公開データ>/data --out models/s2cloud_v1   （s2cloud_train.py の説明）
"""
import os, json
import numpy as np

MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")
DEFAULT = "s2cloud_v1"


def tag(meta, near_override=None):
    """しきい値の決め方の短い名前（方式の名前・置き場に入れる。変えたら作り直しになるように）"""
    t = float(meta.get("threshold", 0.5)); n = float(near_override or meta.get("near_threshold", t))
    out = f"n{n:.2f}" if t == 0.5 else f"t{t:.2f}n{n:.2f}"
    for sv in meta.get("season") or []:
        ms = sv.get("months", [])
        out += f"-m{min(ms)}-{max(ms)}t{float(sv.get('threshold', t)):.2f}n{float(sv.get('near_threshold', n)):.2f}" if ms else ""
    return out
# 計算できる特徴（名前 → 作り方は features）。判定器ごとに、このうち使うものを .json の feats に持つ
KNOWN = ("B02", "B03", "B04", "B08", "B11", "B12", "NDVI", "SCL", "b02m5", "b02s5", "doy")
BANDS10 = ("B02", "B03", "B04", "B08")
BANDS20 = ("B11", "B12")


def doy_of(date):
    """年の中の位置（1月1日 = 0、1年 = 1）"""
    import datetime as dt
    d = dt.date.fromisoformat(date)
    return (d - dt.date(d.year, 1, 1)).days / 365.0


def features(R, date, feats, mask=None):
    """R: {帯の名前: 整数の配列（L2A の値。B11・B12・SCL は 10m に広げたもの）}。戻り値 F[H, W, len(feats)]（float32）。
    mask（bool, H×W）を渡すと、その画素だけの F[N, len(feats)]（速い。値は同じ）。
    学習（probe/s2_open/cs_feat.py）と同じ式: 反射率 = 値/10000、NDVI = (B08−B04)/(B08+B04)、青の 5×5 の平均とばらつき（画面全体で計算してから選ぶ）"""
    from scipy import ndimage
    H, W = R["SCL"].shape
    pick = (lambda a: a) if mask is None else (lambda a: a[mask])
    ref = lambda b: pick(R[b]).astype("float32") / 1e4
    out = []; m5 = None
    for f in feats:
        if f in BANDS10 + BANDS20:
            a = ref(f)
        elif f == "NDVI":
            b4, b8 = ref("B04"), ref("B08"); a = (b8 - b4) / np.maximum(b8 + b4, 1e-6)
        elif f == "SCL":
            a = pick(R["SCL"]).astype("float32")
        elif f in ("b02m5", "b02s5"):
            if m5 is None:
                b2 = R["B02"].astype("float32") / 1e4; m5 = ndimage.uniform_filter(b2, 5); q5 = ndimage.uniform_filter(b2 * b2, 5)
            a = pick(m5) if f == "b02m5" else np.sqrt(np.maximum(pick(q5) - pick(m5) * pick(m5), 0))
        elif f == "doy":
            a = np.full(pick(np.empty((H, W), bool)).shape, doy_of(date), "float32")
        else:
            raise ValueError(f"知らない特徴 {f}")
        out.append(a.astype("float32", copy=False))
    return np.stack(out, -1)


def bands_needed(feats):
    """判定器の特徴を作るのに読む帯（NDVI 用の B04・B08 と SCL は、いつも読む）"""
    need = {"B04", "B08", "SCL"}
    for f in feats:
        if f in BANDS10 + BANDS20:
            need.add(f)
        elif f in ("b02m5", "b02s5"):
            need.add("B02")
    return need


class CloudModel:
    """雲の確からしさ（0〜1）。threshold 以上を雲とする"""
    def __init__(self, name=DEFAULT, model_dir=MODEL_DIR):
        import pickle
        base = name if os.path.isabs(name) else os.path.join(model_dir, name)
        with open(base + ".json", encoding="utf-8") as f:
            self.meta = json.load(f)
        with open(base + ".pkl", "rb") as f:
            self.model = pickle.load(f)
        self.name = os.path.basename(base); self.feats = self.meta["feats"]; self.threshold = float(self.meta.get("threshold", 0.5))
        try:                                   # 学習したときと版が違うと、読めても結果が変わりうる（requirements-open.txt で固定する）
            import sklearn
            if sklearn.__version__ != self.meta.get("sklearn"):
                print(f"注意: 雲の判定器 {self.name} は scikit-learn {self.meta.get('sklearn')} で作った（いまは {sklearn.__version__}）", flush=True)
        except ImportError:
            pass

    def bands(self):
        return bands_needed(self.feats)

    def thresholds(self, date, near_default=None):
        """その日の (画素を使うかのしきい値, 雲の近くの起点のしきい値)。.json の season（月ごとの上書き）があれば、その月の値
        （判定器は 5〜9月だけで学習したので、ほかの月は Earth Engine の判断と食い違いがつり合う値にする。probe/s2_open/open_season.py）"""
        th = self.threshold; nth = near_default if near_default is not None else float(self.meta.get("near_threshold", th))
        m = int(date[5:7])
        for sv in self.meta.get("season") or []:
            if m in sv.get("months", []):
                th = float(sv.get("threshold", th)); nth = float(sv.get("near_threshold", nth))
        return th, nth

    def prob(self, R, date, mask):
        """mask（bool, H×W）の画素の雲の確からしさ。ほかは 0"""
        H, W = mask.shape
        p = np.zeros((H, W), "float32")
        if mask.any():
            p[mask] = self.model.predict_proba(features(R, date, self.feats, mask))[:, 1]
        return p
