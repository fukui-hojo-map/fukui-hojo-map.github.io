"""Earth Engine を使わない Sentinel-2 の計算（AWS の公開データ sentinel-cogs。アカウント・鍵は要らない）。

update_ndvi.py の GEEBackend と同じ口（list_dates・images・stats）と、pixels.py の GEEPixels と同じ口（paddy・ndvi）を持つ。
  python pipeline/update_ndvi.py --backend open ...     （既定は config の s2_backend＝gee。本番はまだ切り替えていない）
  python pipeline/s2open.py prefill --start 2025-05-01 --end 2025-10-01 [--cells ...]   （日ごと・ブロックごとにまとめて先に計算する。全県の作り直し用）
  公開中の ndvi.json との比べは probe/s2_open/open_compare.py（手元の確かめ用）

考え方（Earth Engine 版 update_ndvi.GEEBackend と同じ結果になるように）:
  画像    : Element 84 の sentinel-cogs（Sentinel-2 L2A の COG。BOA の +1000 は適用ずみ）。S3 の一覧（タイル・月ごと）で画像を探し、
            画像ごとの説明（STAC の JSON: 撮影時刻・写っている範囲・処理の版）を読む。同じタイル・同じ撮影（datatake）が2つあれば新しい処理の方
  格子    : UTM（福井は 53帯 EPSG:32653）の 10m。タイルの左上は画像のヘッダから読む（タイルごとに 20m ずれるので、式で決めない。
            決め打ちで 2画素ずれ、7月に −0.03 の誤りになったことがある）。どのタイルも同じ 10m・20m の格子に乗る
  雲      : Cloud Score+ の判断を学習した判定器（s2cloud.py）。田（筆ポリゴン）の画素にだけ使う。SCL 0（データなし）・1（飽和）・11（雪）も除く
  合成    : 同じ日の画像（隣のタイル・別の衛星）は、撮影（datatake）の時刻の順、同じ撮影ならタイルの名前の順に重ね、あとの画像の使える画素を上にする
            （Earth Engine の mosaic は system:index の順＝撮影の時刻・タイルの名前の順で、後ろが上。タイルごとに大気の補正が違い、同じ日・同じ所でも
            NDVI が 0.02〜0.08 違うことがあるので、重ね方を合わせないと値がずれる。試し: 時刻の順では公開値との差 +0.027、名前の順で +0.002）
  雲の近く: 雲の印（SCL 8・9・10、または田の画素で判定器が雲）を 30m の格子（UTM の原点からの格子）に集め、4分の1以上の格子から
            cloud_near_m（300m）以内の画素を使わない（update_ndvi.GEEBackend.day と同じ）
  区画平均: 区画の内側（inner.geojson）が画素にかかる割合（10m 画素を 5×5 に細かくして数える）で重みをつけた平均。数＝重みの合計（10m 画素の数）
処理の単位: 9.6km 四方のブロック（周りに 480m の余白）× 撮影日。COG は 10km 四方のまとまりで圧縮されているので、2km のセルごとに読むと
  同じまとまりを何度も読む（約25倍）。ブロック・日ごとの区画の合計を、置き場（config の s2open_cache、公開しない）に保存し、ほかのセルはそれを使う。
"""
import os, re, sys, json, math, time, threading, hashlib, glob, argparse, urllib.parse, urllib.request, datetime as dt
from concurrent.futures import ThreadPoolExecutor
import numpy as np

S3 = "https://sentinel-cogs.s3.us-west-2.amazonaws.com/"
PREFIX = "sentinel-s2-l2a-cogs"
RES = 10
BLOCK_PX = 960            # 9.6km。10・20・30m のどれでも割り切れる（20m の帯・30m の雲の格子が、ブロックの端でずれない）
HALO_PX = 48              # 480m。雲の近く（300m）＋ 30m の格子 ＋ 5×5 の計算に足りる
WIN_PX = BLOCK_PX + 2 * HALO_PX
K = 5                     # 区画の重み: 10m 画素を 5×5 に細かくして数える
VERSION = "o5"            # 計算のしかたの版。変えると置き場の結果を作り直す（o2: 同じ日の画像の重ね方を Earth Engine に合わせた。o3: 60m の晴れ・雲の数も置く。o4: 60m の田の画素の数・晴れの数も置く。o5: 区画の集め方を Earth Engine と同じに）
# 区画の集め方（2026-11-05、利用者の指示「GEE と同じやり方に合わせて」）: Earth Engine の reduceRegions(scale=10) は、その日の画像（mosaic）の既定の投影の
#   経緯度（EPSG:4326）の格子で集める。1マス = 10m ÷ 赤道の1度の長さ（福井で 東西 8.1m × 南北 10m）。マスの値は、マスの中心にいちばん近い 10m 画素。
#   平均は、マスに区画がかかる割合の重みつき。数（有効画素率に使う）は、区画に少しでもかかるマスの数（重みなし）。
#   公開中の値と比べて確かめた（probe/s2_open/agg_compare.py: 晴れの 12,541 区画・日で 差 0.0000、平均絶対 0.002、95%点 0.005。前の集め方は 0.013・0.056）
EE_DEG = 10 / 111319.49079327357
EE_K = 5                  # マスにかかる割合を、5×5 の点で数える
MIN_FIELD_PX = 300        # セルの晴れの割合を、田の画素だけで数えるのに要る田の画素の数（3ha）
Q = 6                     # 晴れ・雲の数を数える格子（60m = 10m の 6×6）。セルの晴れの割合（生育傾向・px-archive）と、まわりの雲の割合（cloudnear）に使う
MGRS_COL = {1: "ABCDEFGH", 2: "JKLMNPQR", 0: "STUVWXYZ"}
MGRS_ROW = "ABCDEFGHJKLMNPQRSTUV"
LAT_BAND = "CDEFGHJKLMNPQRSTUVWX"
TILE_M = 109800           # タイルの幅（m）
SUPPORTED = ("s2",)       # この backend が、鍵なしで計算できる観測元（Sentinel-2）


def supported():
    """この backend で計算できる観測元。Sentinel-2 は鍵なし（AWS の公開データ）。
    レーダー（Sentinel-1）は Copernicus Data Space の鍵（CDSE_CLIENT_ID・CDSE_CLIENT_SECRET）があるときだけ（s1cdse.py）。
    Landsat はまだ（NASA Earthdata の鍵が要る）。計算できない観測元は、ndvi.json にある分をそのまま残す"""
    return SUPPORTED + (("s1",) if os.environ.get("CDSE_CLIENT_ID") and os.environ.get("CDSE_CLIENT_SECRET") else ())


def read_json(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_json(path, obj):
    """書きかけのファイルを残さない（別の処理が同時に読んでもよいように、書いてから入れ替える）"""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.{os.getpid()}.{threading.get_ident()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, separators=(",", ":"))
    os.replace(tmp, path)


def _http(url, timeout=60, tries=5):
    """S3・説明の JSON を読む（一時的な 5xx・切断は待ってやり直す）"""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            code = getattr(e, "code", None)
            if k == tries - 1 or (code is not None and code < 500 and code != 429):
                raise
            time.sleep(2 * 2 ** k)


def _gdal_env():
    import rasterio
    opt = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MULTIRANGE="YES",
               GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2",
               CPL_VSIL_CURL_CACHE_SIZE=str(int(os.environ.get("S2OPEN_CURL_CACHE_MB", "1500")) << 20))   # 読んだ所の置き場（ファイルを開き直しても使える）。隣のブロックと重なる所を読み直さない
    ca = os.environ.get("GDAL_HTTP_CAINFO") or os.environ.get("CURL_CA_BUNDLE")
    if ca:
        opt["GDAL_HTTP_CAINFO"] = ca
    return rasterio.Env(**opt)


def utm_epsg_of(lon):
    return 32600 + int((lon + 180) // 6) + 1


def mgrs_tile(zone, lat, e, n):
    """UTM の座標が入る 100km 四方の MGRS の名前（北半球）"""
    col = MGRS_COL[zone % 3][int(e // 100000) - 1]
    row = MGRS_ROW[(int(n // 100000) + (0 if zone % 2 else 5)) % 20]
    return f"{zone}{LAT_BAND[int((lat + 80) // 8)]}{col}{row}"


# ---------------- 画像の一覧 ----------------
class Catalog:
    """タイル・月ごとの S3 の一覧と、画像ごとの説明（STAC の JSON）。終わった月の一覧と説明は置き場に残す"""
    def __init__(self, cache_dir):
        self.dir = cache_dir; self.lock = threading.Lock(); self._month = {}; self._item = {}
        os.makedirs(os.path.join(cache_dir, "items"), exist_ok=True); os.makedirs(os.path.join(cache_dir, "list"), exist_ok=True)
        self.tiles_p = os.path.join(cache_dir, "tiles.json"); self._tiles = read_json(self.tiles_p) or {}

    def tile_ul(self, tile):
        """タイルの左上（UTM）。画像のヘッダから読む。タイルがない（画像が一つもない）ときは None"""
        if tile in self._tiles:
            return self._tiles[tile]
        ul = None
        today = dt.date.today()
        for back in (40, 80, 400):
            d = today - dt.timedelta(days=back)
            scs = self.month(tile, d.year, d.month)
            if scs:
                import rasterio
                with _gdal_env(), rasterio.open(f"{scs[0]}/B04.tif") as s:
                    ul = [s.transform.c, s.transform.f]
                break
        with self.lock:
            self._tiles[tile] = ul; write_json(self.tiles_p, self._tiles)
        return ul

    def tiles_for(self, bbox):
        """範囲（経緯度）に写るタイル: 範囲の中を 5km おきに見た 100km 四方と、その西・南の隣（タイルは 9.8km 重なる）"""
        from pyproj import Transformer
        zone = utm_epsg_of((bbox[0] + bbox[2]) / 2) - 32600
        tr = Transformer.from_crs(4326, 32600 + zone, always_xy=True)
        cand = set()
        for lon in np.linspace(bbox[0], bbox[2], max(2, int((bbox[2] - bbox[0]) / 0.05) + 1)):
            for lat in np.linspace(bbox[1], bbox[3], max(2, int((bbox[3] - bbox[1]) / 0.045) + 1)):
                e, n = tr.transform(lon, lat)
                for de, dn in ((0, 0), (-100000, 0), (0, 100000), (-100000, 100000)):
                    cand.add(mgrs_tile(zone, lat, e + de, n + dn))
        (x0, y0), (x1, y1) = tr.transform(bbox[0], bbox[1]), tr.transform(bbox[2], bbox[3])
        xs = [tr.transform(bbox[i], bbox[j])[0] for i in (0, 2) for j in (1, 3)]; ys = [tr.transform(bbox[i], bbox[j])[1] for i in (0, 2) for j in (1, 3)]
        out = []
        for t in sorted(cand):
            ul = self.tile_ul(t)
            if ul and ul[0] < max(xs) and ul[0] + TILE_M > min(xs) and ul[1] - TILE_M < max(ys) and ul[1] > min(ys):
                out.append(t)
        return out

    def month(self, tile, y, m):
        """その月のタイルの画像の場所（URL の前の部分）の一覧"""
        key = (tile, y, m)
        if key in self._month:
            return self._month[key]
        path = os.path.join(self.dir, "list", f"{tile}_{y}_{m:02d}.json")
        done = dt.date(y, m, 1) + dt.timedelta(days=75) < dt.date.today()   # 終わって1か月以上たった月は、もう増えない
        if done and os.path.exists(path):
            out = read_json(path)
        else:
            pre = f"{PREFIX}/{tile[:2]}/{tile[2]}/{tile[3:]}/{y}/{m}/"; out = []; tok = ""
            while True:
                u = S3 + "?list-type=2&delimiter=/&prefix=" + urllib.parse.quote(pre) + (f"&continuation-token={urllib.parse.quote(tok)}" if tok else "")
                x = _http(u).decode()
                out += [S3 + p for p in re.findall(r"<Prefix>(" + re.escape(pre) + r"[^<]*)/</Prefix>", x)]
                mm = re.search(r"<NextContinuationToken>([^<]*)<", x); tok = mm.group(1) if mm else ""
                if not tok:
                    break
            out = sorted(out)
            if done:
                write_json(path, out)
        with self.lock:
            self._month[key] = out
        return out

    def item(self, url):
        """画像の説明: {id, tile, t（ミリ秒）, date, url, geom（写っている範囲 GeoJSON）, baseline, datatake, updated}"""
        sid = url.rstrip("/").split("/")[-1]
        if sid in self._item:
            return self._item[sid]
        path = os.path.join(self.dir, "items", sid + ".json")
        it = read_json(path)
        if not it:
            j = json.loads(_http(f"{url}/{sid}.json"))
            p = j["properties"]; t = dt.datetime.fromisoformat(p["datetime"].replace("Z", "+00:00"))
            it = {"id": sid, "tile": f"{p.get('mgrs:utm_zone')}{p.get('mgrs:latitude_band')}{p.get('mgrs:grid_square')}", "url": url,
                  "t": int(t.timestamp() * 1000), "date": t.strftime("%Y-%m-%d"), "geom": j["geometry"],
                  "baseline": p.get("s2:processing_baseline") or "", "datatake": p.get("s2:datatake_id") or sid, "updated": p.get("updated") or ""}
            write_json(path, it)
        with self.lock:
            self._item[sid] = it
        return it

    def scenes(self, tiles, start, end, bbox=None):
        """start 以上 end 未満の日の画像（bbox と重なるもの）。重ねる順（mosaic_key）。同じ画像の処理し直しは新しい方だけ（dedup）"""
        from shapely.geometry import shape, box
        a, b = dt.date.fromisoformat(start), dt.date.fromisoformat(end) - dt.timedelta(days=1)
        urls = []
        y, m = a.year, a.month
        while (y, m) <= (b.year, b.month):
            for t in tiles:
                for u in self.month(t, y, m):
                    mm = re.search(r"_(\d{8})_\d+_L2A$", u)
                    if mm and start <= f"{mm.group(1)[:4]}-{mm.group(1)[4:6]}-{mm.group(1)[6:]}" < end:
                        urls.append(u)
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
        with ThreadPoolExecutor(8) as ex:
            items = list(ex.map(self.item, urls))
        bb = box(*bbox) if bbox else None
        out = [it for it in dedup(items) if start <= it["date"] < end and (bb is None or shape(it["geom"]).intersects(bb))]
        return sorted(out, key=mosaic_key)


def read_window(cat, it, x0, y1, W, H, bands):
    """画像 it の窓（左上 x0・y1、幅 W・高さ H 画素、10m）を読む（タイルの外は 0）。20m の帯（B11・B12・SCL）は 10m に広げる（最も近い画素）"""
    import rasterio
    from rasterio.windows import Window
    ul = cat.tile_ul(it["tile"])
    c0 = int(round((x0 - ul[0]) / RES)); r0 = int(round((ul[1] - y1) / RES))
    def one(b):
        half = b in ("B11", "B12", "SCL")
        cc, rr, ww, hh = (c0 // 2, r0 // 2, W // 2, H // 2) if half else (c0, r0, W, H)
        for k in range(5):
            try:   # やり直すときは URL の後ろを変える（GDAL は一度失敗したファイルを「ない」と覚えてしまうため。S3 は ?以降を無視する）
                with _gdal_env(), rasterio.open(f"{it['url']}/{b}.tif" + (f"?retry={k}" if k else "")) as s:
                    out = np.zeros((hh, ww), s.dtypes[0])
                    ca, ra = max(cc, 0), max(rr, 0); cb, rb = min(cc + ww, s.width), min(rr + hh, s.height)
                    if cb > ca and rb > ra:
                        out[ra - rr:rb - rr, ca - cc:cb - cc] = s.read(1, window=Window(ca, ra, cb - ca, rb - ra))
                break
            except Exception:
                if k == 4:
                    raise
                time.sleep(3 * 2 ** k)
        return b, (np.repeat(np.repeat(out, 2, 0), 2, 1) if half else out)
    with ThreadPoolExecutor(len(bands)) as ex:
        return dict(ex.map(one, sorted(bands)))


def dedup(items):
    """同じタイル・同じ撮影（datatake）で、写っている範囲もほぼ同じものは「処理し直し」とみなし、新しい処理の方だけ残す。
    範囲が違うもの（撮影の帯の切れ目で、1つのタイルが2つに分かれて届く。例: 2026-03-30 の 53SNV は面積 0.24 と 0.97）は両方残して重ねる
    （片方だけにすると、その日の田の一部が「データなし」になる）"""
    from shapely.geometry import shape
    groups = {}
    for it in items:
        groups.setdefault((it["tile"], it["datatake"]), []).append(it)
    out = []
    for g in groups.values():
        keep = []
        for it in sorted(g, key=lambda x: (x["baseline"], x["updated"]), reverse=True):   # 新しい処理から
            a = shape(it["geom"])
            same = any(a.symmetric_difference(shape(k["geom"])).area <= 0.01 * max(a.union(shape(k["geom"])).area, 1e-12) for k in keep)
            if not same:
                keep.append(it)
        out += keep
    return out


def mosaic_key(it):
    """同じ日の画像の重ね方（あとが上）: 撮影（datatake）の時刻、同じ撮影ならタイルの名前（Earth Engine の system:index の順）"""
    m = re.search(r"(\d{8}T\d{6})", it.get("datatake") or "")
    return (m.group(1) if m else dt.datetime.fromtimestamp(it["t"] / 1000, dt.timezone.utc).strftime("%Y%m%dT%H%M%S"), it["tile"], it.get("t", 0), it.get("id", ""))


def hits(cat, it, x0, y1, W, H):
    """画像のタイルが窓にかかるか"""
    ul = cat.tile_ul(it["tile"])
    return bool(ul) and ul[0] < x0 + W * RES and ul[0] + TILE_M > x0 and ul[1] - TILE_M < y1 and ul[1] > y1 - H * RES


# ---------------- 区画 ----------------
class Parcels:
    """県全体の区画（セルの inner.geojson ＝重みの計算、parcels.geojson ＝判定器を使う田の画素）を UTM で持つ"""
    def __init__(self, data_dir, epsg):
        import shapely
        from shapely.geometry import shape
        from pyproj import Transformer
        tr = Transformer.from_crs(4326, epsg, always_xy=True)
        f = lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1]))
        self.pid = []; self.cell = []; inner = []; full = []; h = hashlib.sha1()
        for cdir in sorted(glob.glob(os.path.join(data_dir, "cells", "*"))):
            cid = os.path.basename(cdir)
            ip, pp = os.path.join(cdir, "inner.geojson"), os.path.join(cdir, "parcels.geojson")
            if not os.path.exists(ip):
                continue
            g = read_json(ip) or {"features": []}
            h.update(cid.encode())
            for ft in g["features"]:
                self.pid.append(ft["properties"]["pid"]); self.cell.append(cid); inner.append(shape(ft["geometry"]))
                h.update(ft["properties"]["pid"].encode()); h.update(json.dumps(ft["geometry"]["coordinates"], separators=(",", ":")).encode())
            for ft in (read_json(pp) or {"features": []})["features"]:
                full.append(shape(ft["geometry"]))
        self.hash = h.hexdigest()[:12]
        self.inner_ll = np.array(inner, dtype=object) if inner else np.array([], dtype=object)   # 経緯度のまま（Earth Engine と同じ集め方の重み）
        self.inner = shapely.transform(self.inner_ll, f) if inner else np.array([], dtype=object)
        self.full = shapely.transform(np.array(full, dtype=object), f) if full else np.array([], dtype=object)
        self.index = {p: i for i, p in enumerate(self.pid)}
        self.tree_i = shapely.STRtree(self.inner); self.tree_f = shapely.STRtree(self.full)
        b = shapely.bounds(self.inner) if len(self.inner) else np.zeros((0, 4))
        self.bounds = b                                       # [x0, y0, x1, y1]（UTM）

    def blocks_of(self, gis):
        """区画（添字）にかかるブロック（bx, by）"""
        out = set()
        for x0, y0, x1, y1 in self.bounds[list(gis)]:
            for bx in range(int(x0 // (BLOCK_PX * RES)), int(x1 // (BLOCK_PX * RES)) + 1):
                for by in range(int(y0 // (BLOCK_PX * RES)), int(y1 // (BLOCK_PX * RES)) + 1):
                    out.add((bx, by))
        return sorted(out)


def block_window(bkey):
    """ブロックの読む窓（余白つき）: (左の x, 上の y, 幅, 高さ)（UTM, m・画素）"""
    bx, by = bkey
    return bx * BLOCK_PX * RES - HALO_PX * RES, (by + 1) * BLOCK_PX * RES + HALO_PX * RES, WIN_PX, WIN_PX


# ---------------- 計算の本体 ----------------
class OpenS2:
    """ブロック・日ごとの計算と、その置き場。update_ndvi の backend（OpenBackend）と圃場内マップ（OpenPixels）が共有する"""
    def __init__(self, cfg, data_dir=None, model=None):
        from s2cloud import CloudModel
        self.cfg = cfg
        self.data_dir = data_dir or os.path.join(cfg["site_dir"], "data")
        self.cache = cfg.get("s2open_cache") or "./state/s2open"
        self.model = CloudModel(model or cfg.get("s2open_model") or "s2cloud_v1")
        self.near = int(cfg.get("cloud_near_m") or 0)
        # 雲の近くの起点にする判定器の雲のしきい値（画素を使うかどうかの threshold とは別。判定器の点々とした誤りで、300m 四方を広く外しすぎないように）
        self._near_override = float(cfg["s2open_near_th"]) if cfg.get("s2open_near_th") else None   # config で雲の近くの起点のしきい値を決め打ちするとき（月ごとの上書きより強い）
        self.near_th = self._near_override or float(self.model.meta.get("near_threshold") or self.model.threshold)
        idx = read_json(os.path.join(self.data_dir, "index.json")) or {}
        self.bbox = idx.get("bbox") or cfg.get("bbox")
        self.epsg = utm_epsg_of((self.bbox[0] + self.bbox[2]) / 2)
        self.cat = Catalog(os.path.join(self.cache, "catalog"))
        self._parcels = None; self._tiles = None; self._lock = threading.Lock(); self._klocks = {}
        self._w = {}; self._field = {}; self._px = {}; self._px_order = []; self._pxmap = {}
        (sm, sd), (em, ed) = [x.split("-") for x in cfg.get("pixel_window", ["05-10", "09-30"])]
        self.px_md = (f"{sm}-{sd}", f"{em}-{ed}")
        self._day = {}                                             # 日 → その日の画像（この実行の間だけ）
        self.stats_n = {"blocks": 0, "scenes": 0, "read_s": 0.0, "calc_s": 0.0, "cache_hit": 0}

    # --- 準備 ---
    @property
    def parcels(self):
        with self._lock:
            if self._parcels is None:
                self._parcels = Parcels(self.data_dir, self.epsg)
            return self._parcels

    @property
    def tiles(self):
        if self._tiles is None:
            self._tiles = self.cat.tiles_for(self.bbox)
        return self._tiles

    def res_dir(self):
        from s2cloud import tag
        return os.path.join(self.cache, "res", f"{VERSION}-{self.model.name}-{tag(self.model.meta, self._near_override)}-near{self.near}-{self.parcels.hash}")

    def scenes_on(self, date):
        if date not in self._day:
            nxt = (dt.date.fromisoformat(date) + dt.timedelta(days=1)).isoformat()
            self._day[date] = self.cat.scenes(self.tiles, date, nxt, self.bbox)
        return self._day[date]

    def prefetch_days(self, start, end):
        """期間の日ごとの画像をまとめて調べておく（一覧の読み込みを月ごと1回に）"""
        by = {}
        for it in self.cat.scenes(self.tiles, start, end, self.bbox):
            by.setdefault(it["date"], []).append(it)
        d = dt.date.fromisoformat(start)
        while d.isoformat() < end:
            self._day.setdefault(d.isoformat(), by.get(d.isoformat(), [])); d += dt.timedelta(days=1)
        return sorted(k for k, v in by.items())

    def _klock(self, key):
        with self._lock:
            return self._klocks.setdefault(key, threading.Lock())

    def weights(self, bkey):
        """ブロックの芯（余白を除く 9.6km 四方）の区画の重み: (区画の添字, 窓の中の画素の番号, 重み, マスの数)。置き場に残す。
        Earth Engine と同じ集め方（上の EE_DEG）: 重み = その画素を「いちばん近い画素」とする経緯度のマスに、区画がかかる割合の和。マスの数 = そのマスの数"""
        if bkey in self._w:
            return self._w[bkey]
        path = os.path.join(self.cache, "weights_ee1", self.parcels.hash, f"{bkey[0]}_{bkey[1]}.npz")
        if os.path.exists(path):
            z = np.load(path); w = (z["g"], z["pix"], z["w"], z["c"])
        else:
            w = self._make_weights(bkey)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = path + f".{os.getpid()}.{threading.get_ident()}.tmp.npz"
            np.savez_compressed(tmp, g=w[0], pix=w[1], w=w[2], c=w[3]); os.replace(tmp, path)
        self._w[bkey] = w
        return w

    def _make_weights(self, bkey):
        """経緯度のマス（EE_DEG）に区画を塗り（EE_K×EE_K の点）、マスの中心に最も近い 10m 画素がこのブロックの芯にあるマスだけを、その画素に集める"""
        import shapely
        from rasterio.features import rasterize
        from rasterio.transform import Affine
        from pyproj import Transformer
        bx, by = bkey; L = BLOCK_PX * RES; x0, y1 = bx * L, (by + 1) * L
        P = self.parcels
        empty = (np.zeros(0, np.int32), np.zeros(0, np.int32), np.zeros(0, np.float32), np.zeros(0, np.int16))
        gi = P.tree_i.query(shapely.box(x0 - 2 * RES, y1 - L - 2 * RES, x0 + L + 2 * RES, y1 + 2 * RES), predicate="intersects")
        if not len(gi):
            return empty
        inv = Transformer.from_crs(self.epsg, 4326, always_xy=True)
        ex = np.linspace(x0, x0 + L, 9); ey = np.linspace(y1 - L, y1, 9)
        X = np.concatenate([ex, ex, np.full(9, x0), np.full(9, x0 + L)]); Y = np.concatenate([np.full(9, y1 - L), np.full(9, y1), ey, ey])
        lo, la = inv.transform(X, Y)
        i0 = math.floor(min(lo) / EE_DEG) - 2; i1 = math.ceil(max(lo) / EE_DEG) + 2; j0 = math.floor(min(la) / EE_DEG) - 2; j1 = math.ceil(max(la) / EE_DEG) + 2
        W, H = (i1 - i0) * EE_K, (j1 - j0) * EE_K
        lab = rasterize(((P.inner_ll[i], int(i) + 1) for i in gi), out_shape=(H, W),
                        transform=Affine(EE_DEG / EE_K, 0, i0 * EE_DEG, 0, -EE_DEG / EE_K, j1 * EE_DEG), fill=0, dtype="int32")
        rr, cc = np.nonzero(lab)
        if not len(rr):
            return empty
        lv = lab[rr, cc].astype(np.int64) - 1; del lab
        nc = (i1 - i0)
        cell = (rr // EE_K).astype(np.int64) * nc + (cc // EE_K)                   # マスの番号
        key, cnt = np.unique(lv * (nc * (j1 - j0)) + cell, return_counts=True)    # （区画, マス）ごとに、中にある点の数
        kg = key // (nc * (j1 - j0)); kc = key % (nc * (j1 - j0))
        ucell, ui = np.unique(kc, return_inverse=True)
        clon = (i0 + (ucell % nc) + 0.5) * EE_DEG; clat = (j1 - (ucell // nc) - 0.5) * EE_DEG   # マスの中心
        ux, uy = Transformer.from_crs(4326, self.epsg, always_xy=True).transform(clon, clat)
        col = np.floor(np.asarray(ux) / RES).astype(np.int64); row = np.floor(np.asarray(uy) / RES).astype(np.int64)   # 中心にいちばん近い 10m 画素（全体の格子。行は北向き）
        inb = (col // BLOCK_PX == bx) & (row // BLOCK_PX == by)                       # その画素がこのブロックの芯にあるマスだけ（ブロックをまたいで二重に数えない）
        pixc = ((BLOCK_PX - 1) - (row - by * BLOCK_PX) + HALO_PX) * WIN_PX + (col - bx * BLOCK_PX + HALO_PX)
        m = inb[ui]
        if not m.any():
            return empty
        g2 = kg[m]; p2 = pixc[ui][m]; w2 = cnt[m] / (EE_K * EE_K)
        k2, inv2 = np.unique(g2 * (WIN_PX * WIN_PX) + p2, return_inverse=True)
        ws = np.bincount(inv2, weights=w2); cs = np.bincount(inv2)
        return ((k2 // (WIN_PX * WIN_PX)).astype(np.int32), (k2 % (WIN_PX * WIN_PX)).astype(np.int32), ws.astype(np.float32), cs.astype(np.int16))

    def field(self, bkey):
        """窓の中の田の画素（parcels.geojson の形を 1画素広げたもの。判定器を使う画素）"""
        if bkey in self._field:
            return self._field[bkey]
        import shapely
        from rasterio.features import rasterize
        from rasterio.transform import Affine
        x0, y1, W, H = block_window(bkey); P = self.parcels
        fi = P.tree_f.query(shapely.box(x0, y1 - H * RES, x0 + W * RES, y1), predicate="intersects")
        from scipy import ndimage
        a = (rasterize(((P.full[i], 1) for i in fi), out_shape=(H, W), transform=Affine(RES, 0, x0, 0, -RES, y1), fill=0, dtype="uint8") > 0
             if len(fi) else np.zeros((H, W), bool))
        a = ndimage.binary_dilation(a)                     # 1画素広げる（圃場内マップの格子とのずれで、田の縁の雲を見落とさない）
        self._field[bkey] = a
        return a

    # --- 1日・1ブロック ---
    def _read(self, it, x0, y1, W, H, bands):
        return read_window(self.cat, it, x0, y1, W, H, bands)

    def _hits(self, it, x0, y1, W, H):
        return hits(self.cat, it, x0, y1, W, H)

    def want_px(self, date):
        """圃場内マップの期間（config の pixel_window。どの年でも）の日か"""
        return self.px_md[0] <= date[5:] <= self.px_md[1]

    def process(self, date, bkey, keep_px=False, keep_nd=False):
        """1日・1ブロックを計算する。戻り値 {"g": 区画の添字, "s": [n, 6]（遠い画素の重みの和, 遠い画素の重み×NDVI の和,
        使える画素の重みの和, 使える画素の重み×NDVI の和, 遠い画素のマスの数, 使える画素のマスの数）, "scenes": 使った画像, "px": 圃場内マップ用の NDVI（芯だけ, uint8。keep_px のとき）}"""
        x0, y1, W, H = block_window(bkey)
        scenes = sorted((it for it in self.scenes_on(date) if self._hits(it, x0, y1, W, H)), key=mosaic_key)   # あとの画像が上
        g, pix, w, wc = self.weights(bkey)
        field = self.field(bkey)
        nd = np.zeros((H, W), np.float32); valid = np.zeros((H, W), bool); cld = np.zeros((H, W), bool)
        seen = np.zeros((H, W), bool); bad = np.zeros((H, W), bool)         # データあり（SCL≠0）・悪い画素（雲・影・巻雲、田は判定器の雲）
        bands = self.model.bands() if field.any() else {"B04", "B08", "SCL"}   # 田のないブロック（まわりの雲の割合だけ）は判定器を使わない
        th, nth = self.model.thresholds(date, self._near_override) if hasattr(self.model, "thresholds") else (self.model.threshold, self.near_th)
        # 上の画像から順に見る（結果は、下から重ねて上を優先するのと同じ）。下の画像は、上で決まらなかった画素だけ判定器を使う（タイルの重なりで速い）
        #   NDVI: いちばん上の「使える」画素の値。雲の印・悪い画素: いちばん上の「データのある」画像の印
        for it in reversed(scenes):
            t0 = time.time(); R = self._read(it, x0, y1, W, H, bands); t1 = time.time()
            scl = R["SCL"]; b4 = R["B04"].astype(np.float32); b8 = R["B08"].astype(np.float32)
            have = (b4 > 0) & (b8 > 0) & (scl != 0)
            ok = have & ~np.isin(scl, (1, 11))
            pr = self.model.prob(R, date, ok & field & ~(valid & seen))   # 判定器の雲の確からしさ（田の画素で、まだ決まっていない所だけ）
            v = ok & ~(pr >= th) & ~valid
            nd = np.where(v, (b8 - b4) / np.maximum(b8 + b4, 1), nd); valid |= v
            new = (scl != 0) & ~seen
            cld = np.where(new, np.isin(scl, (8, 9, 10)) | (pr >= nth), cld)
            bad = np.where(new, np.isin(scl, (3, 8, 9, 10)) | (pr >= th), bad); seen |= new
            self.stats_n["scenes"] += 1; self.stats_n["read_s"] += t1 - t0; self.stats_n["calc_s"] += time.time() - t1
        dat = seen
        far = self.far(cld) if self.near else np.ones((H, W), bool)
        vv = valid.ravel()[pix]; ff = far.ravel()[pix]; nv = nd.ravel()[pix]
        u, inv = np.unique(g, return_inverse=True)
        S = np.zeros((len(u), 6), np.float64)
        wc = wc.astype(np.float64)
        for j, a in enumerate((w * vv * ff, w * vv * ff * nv, w * vv, w * vv * nv, wc * vv * ff, wc * vv)):   # 5・6列目: 使えたマスの数（Earth Engine の count と同じ。重みなし）
            S[:, j] = np.bincount(inv, weights=a, minlength=len(u))
        core = (slice(HALO_PX, HALO_PX + BLOCK_PX),) * 2; n = BLOCK_PX // Q
        cnt = lambda a: a[core].reshape(n, Q, n, Q).sum(axis=(1, 3)).astype(np.uint8)
        r = {"g": u.astype(np.int32), "s": S.astype(np.float32), "scenes": [it["id"] for it in scenes], "q": np.stack([cnt(dat), cnt(bad), cnt(field), cnt(field & valid)])}   # 田の画素・田の晴れ（判定器）の数も
        if keep_nd:   # 調べる用（probe/s2_open/agg_compare.py）: 窓の NDVI（float32）・使える画素・雲から遠い画素
            r["nd"] = nd; r["valid"] = valid; r["far"] = far
        if keep_px:   # Earth Engine の pixels.GEEPixels.ndvi と同じ符号: NDVI −0.2〜1.0 → 1〜255、雲・なし → 0（雲の近くは除かない）
            r["px"] = np.where(valid, np.round((np.clip(nd, -0.2, 1.0) + 0.2) / 1.2 * 254 + 1), 0).astype(np.uint8)[core]
        self.stats_n["blocks"] += 1
        return r

    def _keep(self, date, bkey, a):
        with self._lock:
            if (date, bkey) not in self._px:
                self._px_order.append((date, bkey))
            self._px[(date, bkey)] = a
            while len(self._px_order) > int(self.cfg.get("s2open_px_keep") or 400):
                self._px.pop(self._px_order.pop(0), None)

    def far(self, cld):
        """雲の印（10m）から、雲の近くでない画素（True）。30m の格子で 4分の1以上が雲の格子から cloud_near_m 以上"""
        from scipy import ndimage
        H, W = cld.shape
        c30 = cld.reshape(H // 3, 3, W // 3, 3).mean(axis=(1, 3)) >= 0.25
        if not c30.any():
            return np.ones((H, W), bool)
        return np.kron(ndimage.distance_transform_edt(~c30) * 30 >= self.near, np.ones((3, 3), bool)).astype(bool)

    def get(self, date, bkey):
        """置き場にあればそれ（使った画像が今と同じとき）、なければ計算して置く。圃場内マップの期間の日は、NDVI の画像も置く"""
        path = os.path.join(self.res_dir(), date, f"{bkey[0]}_{bkey[1]}.npz")
        x0, y1, W, H = block_window(bkey)
        ids = [it["id"] for it in sorted(self.scenes_on(date), key=mosaic_key) if self._hits(it, x0, y1, W, H)]
        px = self.want_px(date)
        with self._klock((date, bkey)):
            if os.path.exists(path):
                z = np.load(path)
                if list(z["scenes"]) == ids and (not px or "px" in z):
                    self.stats_n["cache_hit"] += 1
                    return {"g": z["g"], "s": z["s"], "scenes": ids, "path": path}
            r = self.process(date, bkey, px)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = path + f".{os.getpid()}.{threading.get_ident()}.tmp.npz"
            np.savez_compressed(tmp, g=r["g"], s=r["s"], q=r["q"], scenes=np.array(r["scenes"], dtype=str), **({"px": r["px"]} if px else {}))
            os.replace(tmp, path)
            if px:
                self._keep(date, bkey, r["px"])
            r["path"] = path
            return r

    def coarse(self, date, bkey):
        """ブロックの芯の 60m 格子の [データありの 10m 画素の数, 悪い画素の数, 田の画素の数, 田の晴れた画素の数]（各 0〜36、[4, 160, 160]）"""
        return np.load(self.get(date, bkey)["path"])["q"]

    def clear_frac(self, date, bbox):
        """範囲（経緯度の四角）の中の、晴れの割合（データなし・雲・影・巻雲・田の判定器の雲は晴れでない。Earth Engine の trend.clear_dates と同じ考え）。
        範囲に田の画素が MIN_FIELD_PX 以上あれば、田の画素だけで数える（田の晴れ＝学習した判定器）。田の外は SCL しかなく、SCL は雲の周りの晴れた所も
        雲・影（3・8・9・10）にしやすい（2025-06-20 の福井市のセル: SCL では 99.9% が雲、判定器では田の 100% が晴れ。Earth Engine も全部の区画を使った）"""
        from pyproj import Transformer
        tr = Transformer.from_crs(4326, self.epsg, always_xy=True)
        xs, ys = tr.transform([bbox[0], bbox[2], bbox[0], bbox[2]], [bbox[1], bbox[1], bbox[3], bbox[3]])
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys); L = BLOCK_PX * RES; c = Q * RES
        ok = tot = fok = ftot = 0
        for bx in range(int(x0 // L), int(x1 // L) + 1):
            for by in range(int(y0 // L), int(y1 // L) + 1):
                cx = bx * L + (np.arange(BLOCK_PX // Q) + 0.5) * c; cy = (by + 1) * L - (np.arange(BLOCK_PX // Q) + 0.5) * c
                mx = (cx >= x0) & (cx < x1); my = (cy >= y0) & (cy < y1)
                if not (mx.any() and my.any()):
                    continue
                q = self.coarse(date, (bx, by)).astype(np.int64)
                sel = np.ix_(my, mx); tot += my.sum() * mx.sum() * Q * Q; ok += int((q[0][sel] - q[1][sel]).sum())
                if len(q) >= 4:
                    ftot += int(q[2][sel].sum()); fok += int(q[3][sel].sum())
        if ftot >= MIN_FIELD_PX:
            return fok / ftot
        return ok / tot if tot else 0.0

    def near_cloud(self, date, lon, lat, r_m=3000.0):
        """点（経緯度）から半径 r_m の円の中の、悪い画素（雲・影・巻雲。田は判定器の雲）の割合。データのない画素は数えない。データがなければ None
        （cloudnear.py の Earth Engine 版: cs_cdf<0.6 または SCL 3・8・9・10 の割合、60m。田の外は判定器を使えないので SCL だけ）"""
        from pyproj import Transformer
        x, y = Transformer.from_crs(4326, self.epsg, always_xy=True).transform(lon, lat)
        L = BLOCK_PX * RES; c = Q * RES; nb = nd_ = 0
        for bx in range(int((x - r_m) // L), int((x + r_m) // L) + 1):
            for by in range(int((y - r_m) // L), int((y + r_m) // L) + 1):
                cx = bx * L + (np.arange(BLOCK_PX // Q) + 0.5) * c; cy = (by + 1) * L - (np.arange(BLOCK_PX // Q) + 0.5) * c
                m = (cx[None, :] - x) ** 2 + (cy[:, None] - y) ** 2 <= r_m ** 2
                if not m.any():
                    continue
                q = self.coarse(date, (bx, by)).astype(np.int64)
                nd_ += int(q[0][m].sum()); nb += int(q[1][m].sum())
        return nb / nd_ if nd_ else None

    def near_plan(self, lon, lat, r_m=3000.0):
        """near_cloud の下ごしらえ: [(ブロック, 円の中の 60m マスの添字（[160×160] を1列にしたもの）)]。日によらないので、1回作って使い回す（cloudnear.py）"""
        from pyproj import Transformer
        x, y = Transformer.from_crs(4326, self.epsg, always_xy=True).transform(lon, lat)
        L = BLOCK_PX * RES; c = Q * RES; out = []
        for bx in range(int((x - r_m) // L), int((x + r_m) // L) + 1):
            for by in range(int((y - r_m) // L), int((y + r_m) // L) + 1):
                cx = bx * L + (np.arange(BLOCK_PX // Q) + 0.5) * c; cy = (by + 1) * L - (np.arange(BLOCK_PX // Q) + 0.5) * c
                m = (cx[None, :] - x) ** 2 + (cy[:, None] - y) ** 2 <= r_m ** 2
                if m.any():
                    out.append(((bx, by), np.flatnonzero(m)))
        return out

    def near_blocks(self, lon, lat, r_m=3000.0):
        """near_cloud が読むブロック（点から半径 r_m の円にかかる四角）"""
        from pyproj import Transformer
        x, y = Transformer.from_crs(4326, self.epsg, always_xy=True).transform(lon, lat); L = BLOCK_PX * RES
        return {(bx, by) for bx in range(int((x - r_m) // L), int((x + r_m) // L) + 1) for by in range(int((y - r_m) // L), int((y + r_m) // L) + 1)}

    def clear_dates(self, bbox, start, end, min_frac):
        """期間内で、範囲の min_frac 以上が晴れている日"""
        days = sorted({it["date"] for it in self.cat.scenes(self.tiles, start, end, bbox)})
        return [d for d in days if self.clear_frac(d, bbox) >= min_frac]

    def px_core(self, date, bkey):
        """ブロックの芯の、その日の圃場内マップ用の NDVI（uint8）。置き場から読む（期間の外の日は、その場で計算する）"""
        a = self._px.get((date, bkey))
        if a is not None:
            return a
        if self.want_px(date):
            r = self.get(date, bkey)
            z = np.load(r["path"]); a = z["px"]
        else:
            a = self.process(date, bkey, True)["px"]
        self._keep(date, bkey, a)
        return a

    def parcel_sums(self, gis, date):
        """区画（添字の並び）の、その日の合計 [n, 6]（ブロックにまたがる区画は足す）。画像がない日は 0"""
        gis = np.asarray(gis, np.int64); out = np.zeros((len(gis), 6), np.float64)
        order = np.argsort(gis); sg = gis[order]
        for bk in self.parcels.blocks_of(gis):
            r = self.get(date, bk)
            if not len(r["g"]):
                continue
            pos = np.searchsorted(sg, r["g"]); pos = np.clip(pos, 0, len(sg) - 1); hit = sg[pos] == r["g"]
            np.add.at(out, order[pos[hit]], r["s"][hit])
        return out

    def px_layer(self, date, g, cell_bbox):
        """圃場内マップの格子 g（EPSG:3857）の、その日の NDVI（uint8, 0=雲・なし）。最も近い 10m 画素の値（Earth Engine の computePixels と同じ）"""
        from pyproj import Transformer
        key = (g["x0"], g["y1"], g["w"], g["h"])
        if key not in self._pxmap:
            tr = Transformer.from_crs(3857, self.epsg, always_xy=True)
            xs = g["x0"] + (np.arange(g["w"]) + 0.5) * g["px"]; ys = g["y1"] - (np.arange(g["h"]) + 0.5) * g["px"]
            X, Y = np.meshgrid(xs, ys); ux, uy = tr.transform(X.ravel(), Y.ravel())
            col = np.floor(ux / RES).astype(np.int64); row = np.floor(uy / RES).astype(np.int64)     # 全体の 10m 格子の列・行（行は北向き）
            L = BLOCK_PX
            self._pxmap[key] = (col // L, row // L, col % L, (L - 1) - (row % L))
        bxs, bys, cx, cy = self._pxmap[key]
        out = np.zeros(len(bxs), np.uint8)
        for bk in sorted(set(zip(bxs.tolist(), bys.tolist()))):
            a = self.px_core(date, bk)
            m = (bxs == bk[0]) & (bys == bk[1])
            out[m] = a[cy[m], cx[m]]
        return out.reshape(g["h"], g["w"])


# ---------------- update_ndvi・pixels の口 ----------------
class OpenBackend:
    """update_ndvi.GEEBackend と同じ口。Sentinel-2（AWS）と、鍵があるときのレーダー（Copernicus Data Space。s1cdse.OpenS1）"""
    def __init__(self, cfg, data_dir=None, engine=None):
        self.cfg = cfg; self.eng = engine or OpenS2(cfg, data_dir)
        self._s1 = None; self._s1_lock = threading.Lock()

    @property
    def s1(self):
        with self._s1_lock:
            if self._s1 is None:
                from s1cdse import OpenS1
                self._s1 = OpenS1(self.eng)
            return self._s1

    def images(self, bbox, start, end, src="s2"):
        if src not in supported():
            return []
        if src == "s1":
            return self.s1.images(bbox, start, end)
        return [{"id": it["id"], "t": it["t"], "geom": it["geom"]} for it in self.eng.cat.scenes(self.eng.tiles, start, end, bbox)]

    def list_dates(self, bbox, start, end, src="s2"):
        if src not in supported():
            return []
        if src == "s1":
            return self.s1.list_dates(bbox, start, end)
        return sorted({it["date"] for it in self.eng.cat.scenes(self.eng.tiles, start, end, bbox)})

    def stats(self, bbox, inner_fc_geojson, dates, src="s2", near=None):
        """GEEBackend.stats と同じ形: near なら {日: {pid: [遠い画素の平均|None, 数, 全部の平均, 数]}}、そうでなければ {日: {pid: [平均, 数]}}"""
        if src not in supported():
            raise ValueError(f"open backend は {src} をまだ計算できない")
        if src == "s1":
            return self.s1.stats(bbox, inner_fc_geojson, dates)
        if near is None:
            near = src == "s2" and self.eng.near > 0
        P = self.eng.parcels
        pids = [f["properties"]["pid"] for f in inner_fc_geojson["features"] if f["properties"]["pid"] in P.index]
        gis = [P.index[p] for p in pids]
        out = {}
        for d in dates:
            out[d] = {}
            if not gis or not self.eng.scenes_on(d):
                continue
            S = self.eng.parcel_sums(gis, d)
            for p, (wf, wfn, wa, wan, cf, ca) in zip(pids, S):
                if wa <= 0:
                    continue                                   # 使える画素がない（Earth Engine では値なしで返らない）
                if near:   # 値は重みつきの平均、数は使えたマスの数（Earth Engine の mean と count）
                    out[d][p] = [round(wfn / wf, 3) if wf > 0 else None, int(round(cf)), round(wan / wa, 3), int(round(ca))]
                else:
                    out[d][p] = [round(wan / wa, 3), int(round(ca))]
        return out


class OpenPixels:
    """pixels.GEEPixels と同じ口（圃場内マップ）"""
    def __init__(self, backend):
        self.b = backend; self.eng = backend.eng

    def paddy(self, fc, g):
        """田の画素（区画の形を塗ったもの、画素の中心が中）: 0/1"""
        from rasterio.features import rasterize
        from rasterio.transform import Affine
        from shapely.geometry import shape
        import shapely
        from pyproj import Transformer
        tr = Transformer.from_crs(4326, 3857, always_xy=True)
        geoms = shapely.transform(np.array([shape(f["geometry"]) for f in fc["features"]], dtype=object),
                                  lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))
        if not len(geoms):
            return np.zeros((g["h"], g["w"]), np.uint8)
        return rasterize(((x, 1) for x in geoms), out_shape=(g["h"], g["w"]), transform=Affine(g["px"], 0, g["x0"], 0, -g["px"], g["y1"]),
                         fill=0, dtype="uint8")

    def ndvi(self, bbox, dates, g):
        return {d: self.eng.px_layer(d, g, bbox) for d in dates if self.eng.scenes_on(d)}


def mask_id(cfg):
    """open backend の雲判定の方式の名前（update_ndvi.mask_id が使う）: 判定器の名前と、雲の近くの起点のしきい値（変えたら作り直しになるように）"""
    from s2cloud import tag
    name = cfg.get("s2open_model") or "s2cloud_v1"
    meta = read_json(os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", name + ".json")) or {}
    return f"s2open-{name}-{tag(meta, float(cfg['s2open_near_th']) if cfg.get('s2open_near_th') else None)}"


# ---------------- 先に計算する（全県の作り直し用） ----------------
def prefill(eng, start, end, cells=None, workers=4, deadline=None, log=print, info=None, extra=None):
    """日ごと・ブロックごとにまとめて計算して置き場に入れる（同じ日の隣のブロックは、読んだ COG のまとまりを使い回せる）。計算したブロック・日の数を返す
    cells: セルの集まり（None = 全部）。{セル: そのセルで計算が要る最初の日} で渡すと、その日より前の日は、そのセルのブロックを計算しない（途中まで済んだ作り直しの続き）
    extra: {セル: 足すブロックの集まり}（区画のないブロック。近くの雲の割合 cloudnear が使う、セルのまわりのブロックなど）
    info: 渡した dict に {"last": 最後まで計算した日（なければ None）, "finished": 期間の終わりまで済んだか} を入れる"""
    P = eng.parcels
    need = cells if isinstance(cells, dict) else None
    by_cell = {}
    for i, c in enumerate(P.cell):
        if cells is None or c in cells:
            by_cell.setdefault(c, []).append(i)
    cell_blocks = {c: set(P.blocks_of(v)) | set((extra or {}).get(c, ())) for c, v in by_cell.items()}
    for c, bl in (extra or {}).items():                  # 区画のないセル
        if c not in cell_blocks and (cells is None or c in cells):
            cell_blocks[c] = set(bl)
    blocks_all = set().union(*cell_blocks.values()) if cell_blocks else set()
    days = eng.prefetch_days(start, end)
    log(f"ブロック {len(blocks_all)}、画像のある日 {len(days)}（{start}〜{end}）、タイル {','.join(eng.tiles)}")
    n = 0; fails = []; last = None; finished = True
    def one(d, bk):
        try:
            eng.get(d, bk); return None
        except Exception as e:              # 1つのブロックの失敗で全体を止めない（置き場に入らないので、使うときにもう一度計算する）
            return (d, bk, str(e)[:200])
    for d in days:
        if deadline and time.time() > deadline:
            log("時間切れ"); finished = False; break
        blocks = blocks_all if need is None else set().union(*[bl for c, bl in cell_blocks.items() if need.get(c, start) <= d])
        if blocks:
            t0 = time.time(); s0 = dict(eng.stats_n)
            with ThreadPoolExecutor(workers) as ex:
                bad = [x for x in ex.map(lambda bk: one(d, bk), sorted(blocks, key=lambda b: (-b[1], b[0]))) if x]   # 北から順（隣を続けて読む）
            n += len(blocks) - len(bad); fails += bad
            log(f"{d}: {time.time() - t0:.0f}秒（ブロック {len(blocks)}・画像 {eng.stats_n['scenes'] - s0['scenes']}・置き場から {eng.stats_n['cache_hit'] - s0['cache_hit']}"
                + (f"・失敗 {len(bad)}: {bad[0][2]}" if bad else "") + "）")
        last = d
    if fails:
        log(f"失敗したブロック・日 {len(fails)}（使うときにもう一度計算する）")
    if info is not None:
        info.update(last=last, finished=finished)
    return n


def main():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from common import load_config
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["prefill", "tiles"])
    ap.add_argument("--start"); ap.add_argument("--end")
    ap.add_argument("--cells", help="カンマ区切りのセル ID（その区画だけ）")
    ap.add_argument("--data", help="公開データの data（既定は config の site_dir/data）")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-minutes", type=float, default=0)
    a = ap.parse_args()
    cfg = load_config()
    eng = OpenS2(cfg, a.data)
    if a.cmd == "tiles":
        print(eng.tiles); return
    t0 = time.time()
    prefill(eng, a.start, a.end, set(a.cells.split(",")) if a.cells else None, a.workers, t0 + a.max_minutes * 60 if a.max_minutes else None,
            log=lambda m: print(time.strftime("%H:%M:%S"), m, flush=True))
    print(f"完了 {(time.time() - t0) / 60:.1f}分", eng.stats_n, flush=True)


if __name__ == "__main__":
    main()
