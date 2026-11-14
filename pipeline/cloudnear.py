"""「近くの雲の割合」を、Sentinel-2 の撮影日ごと・セルごとに出して cloudnear.json に保存する。
画面は、この値が 99% 以上の日の Sentinel-2 の観測を使わない（雲の切れ間だけが晴れて見えている日は、うすい雲・巻雲の影響が残りやすいため。
論文の解析 T18 で、何もしない状態から、観測の 89.3%（6〜9月）を残したまま日ごとの共通の揺れ（SD_D）を 21% 減らせた規則。画面の自動除外と併用すると、観測 86.1% で SD_D −29%）。NDVI の値は一切使わず、値の補正もしない。

  python pipeline/cloudnear.py --max-minutes 120 --workers 4

近くの雲の割合（保存は百分率の切り捨て。割合が 0.99 以上 ⇔ 保存値が 99 以上）: セルの中心から半径 3km の円の中で、Sentinel-2 の画素（60m に集計）のうち、次のどれかに当てはまる画素の割合。
  Cloud Score+ の cs_cdf が 0.6 未満（NDVI の計算で雲として除く基準と同じ）、または SCL が 3（雲の影）・8（雲・中確率）・9（雲・高確率）・10（巻雲）。
  SCL が 0（データなし）の画素は数えない。その日の画像がない円は値なし（null）。値なしの日は、画面では使う側に倒す（観測を捨てない）。

cloudnear.json の形式（v1）:
  {"v": 1, "mask": "bad2-cs0.60-r3", "r": 3, "cells": ["<セルID>", ...],
   "d": {"2025-06-03": [<セルごとの割合(%)の整数 or null>, ...], ...}}     # 並びは cells と同じ
  日付は UTC の撮影日（ndvi.json の dates と同じ）。まだ計算していない日は d にない。
画面が使うのは 6/1〜9/30 の日だけ（評価した時期。site/index.html の CN_FROM/CN_TO）。ほかの時期も、先に計算して保存しておく。
差分更新: 計算済みの日は飛ばす。直近 REDO_DAYS 日は、隣のタイルの画像があとから届くことがあるので、毎回計算し直す。
セルの並びが変わったら、全部計算し直す。時間切れ・エラーでも、計算できた日までを保存する。
方式（mask）が変わったとき（2026-11 の Earth Engine → open など）は、cloudnear_next.json に全部の日を作り、そろってから入れ替える。
そのとき、画面のしきい値を前の方式に合わせ直して "th" に書く（calibrate。画面は th があればそれ、なければ 99）。
--backend は、付けなければ config の s2_backend（なければ open ＝ Earth Engine を使わない。田は学習した判定器、田の外は SCL）。
"""
import argparse, os, sys, json, time, datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed
from common import load_config, write_json, read_json, with_retry, pick_backend

RADIUS_KM = 3.0
SCALE_M = 60
CS_MIN = 0.6
MASK = f"bad2-cs{CS_MIN:.2f}-r{RADIUS_KM:g}"
REDO_DAYS = 10
SAVE_EVERY = 8                # この日数ごとに保存（途中で落ちても進捗を失わない）
TH = 99                       # 画面が「ほぼ雲」とみなす下限（%）のもと。方式を変えたときは、前の方式と印の数がそろう値に合わせ直して "th" に書く（calibrate）
SEASON = ("06-01", "09-30")   # 画面がこの決まりを使う時期（site/index.html の CN_FROM・CN_TO）


def calibrate(old, new, ids):
    """方式を変えたとき、画面のしきい値を合わせ直す: 6〜9月の、両方に値があるセル・日で、「前の方式で TH 以上」の数に、新しい方式の「th 以上」の数がいちばん近い th。
    戻り値 (th, {"from": 前の方式, "n": 比べたセル・日, "old": 前の印, "new": 新しい印, "both": 両方})。比べるものがなければ (None, None)"""
    if not old or old.get("v") != 1 or not old.get("d"):
        return None, None
    oi = {c: k for k, c in enumerate(old.get("cells", []))}
    pairs = []
    for d, v in new.items():
        if not (SEASON[0] <= d[5:] <= SEASON[1]) or d not in old["d"]:
            continue
        o = old["d"][d]
        for k, c in enumerate(ids):
            j = oi.get(c)
            if j is None or j >= len(o) or o[j] is None or v[k] is None:
                continue
            pairs.append((o[j], v[k]))
    n_old = sum(1 for a, _ in pairs if a >= TH)
    if len(pairs) < 1000 or not n_old:
        return None, None
    best = min(range(80, 101), key=lambda t: (abs(sum(1 for _, b in pairs if b >= t) - n_old), abs(t - TH)))
    return best, {"from": old.get("mask"), "n": len(pairs), "old": n_old, "new": sum(1 for _, b in pairs if b >= best),
                  "both": sum(1 for a, b in pairs if a >= TH and b >= best)}


def init_ee():
    import ee
    sa, key, project = os.environ.get("EE_SERVICE_ACCOUNT"), os.environ.get("EE_PRIVATE_KEY"), os.environ.get("EE_PROJECT")
    if sa and key:
        ee.Initialize(ee.ServiceAccountCredentials(sa, key_data=key), project=project)
    else:                                         # ローカルで earthengine authenticate 済みの場合
        ee.Initialize(project=project)
    return ee


def list_dates(ee, region, start, end):
    """期間内の撮影日（UTC）。Cloud Score+ が出ている画像がある日だけ（NDVI の計算と同じ条件）"""
    csp = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED").filterBounds(region).filterDate(start, end)
    ids = csp.aggregate_array("system:index")
    col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(region).filterDate(start, end).filter(ee.Filter.inList("system:index", ids))
    ts = with_retry(lambda: col.aggregate_array("system:time_start").getInfo())
    return sorted({dt.datetime.fromtimestamp(t / 1000, dt.timezone.utc).strftime("%Y-%m-%d") for t in ts})


def day_image(ee, region, d):
    """その日の「悪い画素」(0/1)。同じ日の画像（隣のタイルなど）は1枚に合成する"""
    d0 = ee.Date(d); d1 = d0.advance(1, "day")
    sr = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate(d0, d1).filterBounds(region).mosaic()
    csp = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED").filterDate(d0, d1).filterBounds(region).mosaic()
    scl, cdf = sr.select("SCL"), csp.select("cs_cdf")
    bad = cdf.lt(CS_MIN).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10))
    return bad.updateMask(scl.neq(0)).rename("bad")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", help="data フォルダ（省略時は config の site_dir/data）")
    ap.add_argument("--max-minutes", type=float, default=120)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--start", help="さかのぼる開始日（省略時は config の history_months）")
    ap.add_argument("--smoke", action="store_true", help="動作確認: 直近の2日だけ（保存しない）")
    ap.add_argument("--backend", choices=["gee", "open"], default=None, help="既定は config の s2_backend（なければ open ＝ Earth Engine を使わない。s2open.py。田の外は SCL だけ）")
    a = ap.parse_args()
    cfg = load_config()
    a.backend = pick_backend(cfg, a.backend)
    data_dir = a.data or os.path.join(cfg["site_dir"], "data")
    index = read_json(os.path.join(data_dir, "index.json"))
    if not index:
        sys.exit("index.json がありません")
    cells = [(c["id"], (c["bbox"][0] + c["bbox"][2]) / 2, (c["bbox"][1] + c["bbox"][3]) / 2) for c in index["cells"]]
    ids = [c[0] for c in cells]
    t0 = time.time()
    def log(m): print(time.strftime("%H:%M:%S"), m, flush=True)

    lons = [c[1] for c in cells]; lats = [c[2] for c in cells]
    rbox = [min(lons) - 0.05, min(lats) - 0.05, max(lons) + 0.05, max(lats) + 0.05]
    today = dt.datetime.now(dt.timezone.utc).date()
    start = a.start or (today - dt.timedelta(days=30 * cfg["history_months"])).isoformat()
    end = (today + dt.timedelta(days=1)).isoformat()
    if a.backend == "open":
        from s2open import OpenS2
        eng = OpenS2(cfg, data_dir); mask = f"open-{eng.model.name}-bad2-r{RADIUS_KM:g}"
        dates = sorted({it["date"] for it in eng.cat.scenes(eng.tiles, start, end, rbox)})
    else:
        ee = init_ee(); mask = MASK
        region = ee.Geometry.Rectangle(rbox)
        dates = with_retry(lambda: list_dates(ee, region, start, end))
    log(f"セル {len(ids)} / 撮影日 {len(dates)}（{start}〜）")

    # 方式が変わったとき（Earth Engine → open など）は、別のファイル cloudnear_next.json に全部の日を作り、そろってから入れ替える
    # （それまでは、画面は前の方式の cloudnear.json を使い続ける。途中の、日が欠けたファイルを公開しない）
    pub = os.path.join(data_dir, "cloudnear.json"); nxt = os.path.join(data_dir, "cloudnear_next.json")
    cur = read_json(pub)
    same = lambda o: bool(o) and o.get("v") == 1 and o.get("mask") == mask and o.get("cells") == ids
    swap = bool(cur) and cur.get("v") == 1 and cur.get("mask") != mask and cur.get("cells") == ids and not a.smoke
    path = nxt if swap else pub
    old = read_json(path) if swap else cur
    keep = {}; extra = {}
    if same(old):
        keep = dict(old["d"])                   # 古い日は消さない（ndvi.json も古い日を残すので、選別の対象をそろえる）
        extra = {k: old[k] for k in ("th", "cal") if k in old}
    elif old and not swap:
        log("セルの並びが変わったので、全部計算し直します")
    if swap:
        log(f"方式が変わったので（{cur.get('mask')} → {mask}）、cloudnear_next.json に作ってから入れ替えます")
        dates = sorted(set(dates) | {d for d in cur["d"] if d >= start})   # 前の方式にある日は、新しい方式でも計算を試みる（画像がなければ値なし）
    elif os.path.exists(nxt) and not a.smoke:
        os.remove(nxt)                          # 設定を戻したときなどに残った、途中の計算は捨てる
    if cfg.get("cloud_near_th"):                # しきい値を config で決めるとき
        extra["th"] = int(cfg["cloud_near_th"])
    redo_from = (today - dt.timedelta(days=REDO_DAYS)).isoformat()
    todo = [d for d in dates if d not in keep or d >= redo_from]
    todo.sort(reverse=True)                     # 新しい日から（時間切れでも、よく見る最近の分は先にそろう）
    if a.backend == "open":                     # 置き場に計算済みの日（この実行の中で update_ndvi が作ったもの）を先に（読むだけで済む。全部の日を作り直すときに、時間切れで進まないことを防ぐ）
        rd = eng.res_dir(); cached = lambda d: os.path.isdir(os.path.join(rd, d))
        todo.sort(key=lambda d: not cached(d))
        if swap and any(os.path.exists(os.path.join(data_dir, "cells", c, "ndvi_next.json")) for c in ids):
            n0 = len(todo); todo = [d for d in todo if cached(d)]   # NDVI の作り直しがまだ途中: 置き場にない日は、次回に update_ndvi がまとめて計算するので、いまは待つ
            log(f"NDVI の計算し直しが途中なので、置き場にある日だけ計算します（{len(todo)}／{n0} 日）")
    if a.smoke: todo = todo[:2]
    log(f"計算する日 {len(todo)}（計算済み {len(keep)}）")
    def finish():
        """入れ替え（全部の日がそろったときだけ）: 画面のしきい値を前の方式に合わせ直して、cloudnear.json に置く"""
        if not swap or any(d not in keep for d in dates):
            return False
        th, cal = calibrate(cur, keep, ids)
        if th is not None and "th" not in extra:
            extra.update({"th": th, "cal": cal})
            log(f"画面のしきい値: {th}%（前の方式の 99% 以上 {cal['old']} セル・日 ↔ 新しい方式 {cal['new']}、両方 {cal['both']}。比べた {cal['n']}）")
        write_json(pub, {"v": 1, "mask": mask, "r": RADIUS_KM, "cells": ids, **extra, "d": dict(sorted(keep.items()))})
        if os.path.exists(nxt):
            os.remove(nxt)
        log("入れ替えました（cloudnear.json は新しい方式）")
        return True
    if not todo:
        if not finish():
            log("新しく計算する日はありません")
        return

    if a.backend == "open":
        plans = [eng.near_plan(lo, la, RADIUS_KM * 1000) for cid, lo, la in cells]   # セルごとの [(ブロック, 円の中のマス)]（日によらない）
        def work(d):
            """その日の、セルごとの割合。ブロックは1回だけ読んで、かかるセルで使い回す（eng.near_cloud と同じ計算）"""
            memo = {}; v = []
            for pl in plans:
                nb = nd_ = 0
                for bk, ix in pl:
                    if bk not in memo:
                        q = eng.coarse(d, bk).astype("int64")
                        memo[bk] = (q[0].ravel(), q[1].ravel())
                    q0, q1 = memo[bk]
                    nd_ += int(q0[ix].sum()); nb += int(q1[ix].sum())
                v.append(nb / nd_ if nd_ else None)
            return d, [None if x is None else int(100 * x + 1e-6) for x in v]
    else:
        fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point(lo, la).buffer(RADIUS_KM * 1000), {"id": cid}) for cid, lo, la in cells])
    def work_gee(d):
        img = day_image(ee, region, d)
        res = with_retry(lambda: img.reduceRegions(fc, ee.Reducer.mean(), SCALE_M, tileScale=4).getInfo()["features"])
        by = {f["properties"]["id"]: f["properties"].get("mean") for f in res}
        return d, [None if by.get(c) is None else int(100 * by[c] + 1e-6) for c in ids]

    def save():
        if a.smoke: return
        write_json(path, {"v": 1, "mask": mask, "r": RADIUS_KM, "cells": ids, **extra, "d": dict(sorted(keep.items()))})

    deadline = t0 + a.max_minutes * 60; done = fail = 0
    pending = iter(todo)
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {}
        def submit():
            if time.time() > deadline: return False
            d = next(pending, None)
            if d is None: return False
            futs[ex.submit(work if a.backend == "open" else work_gee, d)] = d; return True
        for _ in range(a.workers): submit()
        while futs:
            for f in as_completed(list(futs)):
                d = futs.pop(f)
                try:
                    _, v = f.result(); keep[d] = v; done += 1
                    log(f"{d}: 計算済み（{done}/{len(todo)}）。最大 {max((x for x in v if x is not None), default='-')}%")
                except Exception as e:
                    fail += 1; log(f"失敗 {d}: {str(e)[:200]}")
                if done and done % SAVE_EVERY == 0: save()
                submit(); break
    save()
    finish()
    left = len(todo) - done - fail
    log(f"完了: 計算 {done} 日・失敗 {fail} 日" + (f"・時間切れで未計算 {left} 日（次回に続き）" if left > 0 else ""))
    if done == 0: sys.exit(1)


if __name__ == "__main__":
    main()
