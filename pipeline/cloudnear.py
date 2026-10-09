"""「近くの雲の割合」を、Sentinel-2 の撮影日ごと・セルごとに出して cloudnear.json に保存する。
画面は、この値が 95% 以上の日の Sentinel-2 の観測を使わない（雲の切れ間だけが晴れて見えている日は、うすい雲・巻雲の影響が残りやすいため。
論文の解析 T18 で、観測の 88.9% を残したまま、日ごとの共通の揺れ（SD_D）を 23% 減らせることを確かめた規則）。NDVI の値は一切使わず、値の補正もしない。

  python pipeline/cloudnear.py --max-minutes 120 --workers 4

近くの雲の割合（保存は百分率の切り捨て。95 以上 ⇔ 割合が 0.95 以上）: セルの中心から半径 3km の円の中で、Sentinel-2 の画素（60m に集計）のうち、次のどれかに当てはまる画素の割合。
  Cloud Score+ の cs_cdf が 0.6 未満（NDVI の計算で雲として除く基準と同じ）、または SCL が 3（雲の影）・8（雲・中確率）・9（雲・高確率）・10（巻雲）。
  SCL が 0（データなし）の画素は数えない。その日の画像がない円は値なし（null）。値なしの日は、画面では使う側に倒す（観測を捨てない）。

cloudnear.json の形式（v1）:
  {"v": 1, "mask": "bad2-cs0.60-r3", "r": 3, "cells": ["<セルID>", ...],
   "d": {"2025-06-03": [<セルごとの割合(%)の整数 or null>, ...], ...}}     # 並びは cells と同じ
  日付は UTC の撮影日（ndvi.json の dates と同じ）。まだ計算していない日は d にない。
画面が使うのは 6/1〜9/30 の日だけ（評価した時期。site/index.html の CN_FROM/CN_TO）。ほかの時期も、先に計算して保存しておく。
差分更新: 計算済みの日は飛ばす。直近 REDO_DAYS 日は、隣のタイルの画像があとから届くことがあるので、毎回計算し直す。
方式（mask）かセルの並びが変わったら、全部計算し直す。時間切れ・エラーでも、計算できた日までを保存する。
"""
import argparse, os, sys, json, time, datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed
from common import load_config, write_json, read_json, with_retry

RADIUS_KM = 3.0
SCALE_M = 60
CS_MIN = 0.6
MASK = f"bad2-cs{CS_MIN:.2f}-r{RADIUS_KM:g}"
REDO_DAYS = 10
SAVE_EVERY = 8                # この日数ごとに保存（途中で落ちても進捗を失わない）


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
    a = ap.parse_args()
    cfg = load_config()
    data_dir = a.data or os.path.join(cfg["site_dir"], "data")
    index = read_json(os.path.join(data_dir, "index.json"))
    if not index:
        sys.exit("index.json がありません")
    cells = [(c["id"], (c["bbox"][0] + c["bbox"][2]) / 2, (c["bbox"][1] + c["bbox"][3]) / 2) for c in index["cells"]]
    ids = [c[0] for c in cells]
    t0 = time.time()
    def log(m): print(time.strftime("%H:%M:%S"), m, flush=True)

    ee = init_ee()
    lons = [c[1] for c in cells]; lats = [c[2] for c in cells]
    region = ee.Geometry.Rectangle([min(lons) - 0.05, min(lats) - 0.05, max(lons) + 0.05, max(lats) + 0.05])
    today = dt.datetime.now(dt.timezone.utc).date()
    start = a.start or (today - dt.timedelta(days=30 * cfg["history_months"])).isoformat()
    end = (today + dt.timedelta(days=1)).isoformat()
    dates = with_retry(lambda: list_dates(ee, region, start, end))
    log(f"セル {len(ids)} / 撮影日 {len(dates)}（{start}〜）")

    path = os.path.join(data_dir, "cloudnear.json")
    old = read_json(path)
    keep = {}
    if old and old.get("v") == 1 and old.get("mask") == MASK and old.get("cells") == ids:
        keep = dict(old["d"])                   # 古い日は消さない（ndvi.json も古い日を残すので、選別の対象をそろえる）
    elif old:
        log("方式かセルの並びが変わったので、全部計算し直します")
    redo_from = (today - dt.timedelta(days=REDO_DAYS)).isoformat()
    todo = [d for d in dates if d not in keep or d >= redo_from]
    todo.sort(reverse=True)                     # 新しい日から（時間切れでも、よく見る最近の分は先にそろう）
    if a.smoke: todo = todo[:2]
    log(f"計算する日 {len(todo)}（計算済み {len(keep)}）")
    if not todo:
        log("新しく計算する日はありません"); return

    fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point(lo, la).buffer(RADIUS_KM * 1000), {"id": cid}) for cid, lo, la in cells])
    def work(d):
        img = day_image(ee, region, d)
        res = with_retry(lambda: img.reduceRegions(fc, ee.Reducer.mean(), SCALE_M, tileScale=4).getInfo()["features"])
        by = {f["properties"]["id"]: f["properties"].get("mean") for f in res}
        return d, [None if by.get(c) is None else int(100 * by[c] + 1e-6) for c in ids]

    def save():
        if a.smoke: return
        write_json(path, {"v": 1, "mask": MASK, "r": RADIUS_KM, "cells": ids, "d": dict(sorted(keep.items()))})

    deadline = t0 + a.max_minutes * 60; done = fail = 0
    pending = iter(todo)
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {}
        def submit():
            if time.time() > deadline: return False
            d = next(pending, None)
            if d is None: return False
            futs[ex.submit(work, d)] = d; return True
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
    left = len(todo) - done - fail
    log(f"完了: 計算 {done} 日・失敗 {fail} 日" + (f"・時間切れで未計算 {left} 日（次回に続き）" if left > 0 else ""))
    if done == 0: sys.exit(1)


if __name__ == "__main__":
    main()
