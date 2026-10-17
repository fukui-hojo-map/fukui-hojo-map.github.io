"""圃場内マップの画素データ（px_<年>.png / px_<年>.json）のうち、公開サイトに置かない古い年を、別のブランチ（px-archive）に作って置く。

公開サイト（GitHub Pages の約1GB の上限に近い）に置く画素データは、今年と昨年だけ（pipeline/pixels.py。それより古い年のファイルは消す）。
そこで、古い年（trend_from 年〜）は px-archive ブランチに置き、画面が raw.githubusercontent.com から読む（site/index.html の pxLoad）。
形式は pixels.py と同じ（cells/<id>/px_<年>.png と .json）。

  python pipeline/px_archive.py --data pages/data --out arc --workers 3 --deadline <UNIX秒>
セルごと・年ごとに、次の順で決める（確定した年だけ。期間が終わって20日たった年）。
  1) 置き場（--out）に、同じ方式・同じ格子の確定ずみのものがあれば、何もしない
  2) 公開中のデータ（--data）に確定ずみの px_<年> があれば、そのまま写す（Earth Engine を使わない。公開サイトから外れる前に取っておく）
  3) どちらにもなければ、Earth Engine で作る（Sentinel-2 の晴れた日。pixels.py と同じ作り方）
1セル・1年ずつ保存するので、時間切れで止まっても次回に続きから作る。置き場からファイルを消すことはない。
--backend fake は、Earth Engine を使わない作り物のデータ（動作確認用）。
"""
import os, shutil, time, threading, datetime as dt
from common import write_json, read_json
from pixels import grid_for, read_layers, write_png, cached_paddy, MIN_CLEAR


def year_window(cfg, y):
    (sm, sd), (em, ed) = [map(int, x.split("-")) for x in cfg.get("pixel_window", ["05-10", "09-30"])]
    return dt.date(y, sm, sd), dt.date(y, em, ed)


def closed_years(cfg, today):
    """確定した年（期間が終わって20日たった年）。trend_from 年から"""
    return [y for y in range(int(cfg.get("trend_from", 2022)), today.year + 1)
            if today > year_window(cfg, y)[1] + dt.timedelta(days=20)]


def parse_years(text, cfg, today):
    if not text:
        return closed_years(cfg, today)
    out = []
    for part in text.split(","):
        a, _, b = part.strip().partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return sorted(set(out))


def is_done(meta, mask, g, win):
    return bool(meta) and bool(meta.get("done")) and meta.get("mask") == mask and meta.get("window") == win \
        and all(meta.get(k) == g[k] for k in g)


def arc_done(cdir, y, mask, g, win):
    """置き場にその年が確定ずみであるか（画像が1日でもあるなら、PNG もあること）"""
    meta = read_json(os.path.join(cdir, f"px_{y}.json"))
    return is_done(meta, mask, g, win) and (not meta.get("dates") or os.path.exists(os.path.join(cdir, f"px_{y}.png")))


class GEEArc:
    """Earth Engine から、晴れた日と画素を取る"""
    def __init__(self, backend):
        from pixels import GEEPixels
        from trend import GEETrend
        self.px = GEEPixels(backend); self.tr = GEETrend(backend, self.px)

    def clear_dates(self, bbox, start, end, min_frac):
        return self.tr.clear_dates(bbox, start, end, min_frac)

    def ndvi(self, bbox, dates, g):
        return self.px.ndvi(bbox, dates, g)

    def paddy(self, fc, g):
        return self.px.paddy(fc, g)


class FakeArc:
    """動作確認用（Earth Engine を使わない）"""
    def __init__(self):
        from trend import FakeTrend
        self.tr = FakeTrend()

    def clear_dates(self, bbox, start, end, min_frac):
        return self.tr.clear_dates(bbox, start, end, min_frac)

    def ndvi(self, bbox, dates, g):
        return self.tr.ndvi(bbox, dates, g)

    def paddy(self, fc, g):
        from trend import labels
        return (labels(fc, g, edges=(0,)) > 0).astype("uint8")


def put_year(arc_cdir, y, win, mask, g, dates, layers):
    """1年分を置き場に書く（確定ずみとして）"""
    png_p, meta_p = os.path.join(arc_cdir, f"px_{y}.png"), os.path.join(arc_cdir, f"px_{y}.json")
    kept = sorted(layers)
    os.makedirs(arc_cdir, exist_ok=True)
    if kept:
        write_png(png_p, [layers[d] for d in kept])
    elif os.path.exists(png_p):
        os.remove(png_p)
    write_json(meta_p, {"v": 1, "mask": mask, "dates": kept, "src_dates": sorted(dates), **g, "window": win, "done": True})
    return len(kept)


def update_cell(get_src, cfg, data_dir, out_dir, cell, years, mask, today, log, deadline=None):
    """1セルの、置き場にない年を足す。足した年数を返す"""
    cid = cell["id"]
    cdir, arc_cdir = os.path.join(data_dir, "cells", cid), os.path.join(out_dir, "cells", cid)
    parcels = read_json(os.path.join(cdir, "parcels.geojson"))
    if not parcels or not parcels.get("features"):
        return 0
    fc = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {}, "geometry": f["geometry"]} for f in parcels["features"]]}
    g = grid_for(cell["bbox"]); made = 0
    for y in years:
        ws, we = year_window(cfg, y); win = [ws.isoformat(), we.isoformat()]
        if arc_done(arc_cdir, y, mask, g, win):
            continue
        if deadline and time.time() > deadline:
            break
        site = read_json(os.path.join(cdir, f"px_{y}.json"))
        if is_done(site, mask, g, win):                          # 公開中のものを写す
            os.makedirs(arc_cdir, exist_ok=True)
            for ext in ("png", "json"):
                p = os.path.join(cdir, f"px_{y}.{ext}")
                if os.path.exists(p):
                    shutil.copyfile(p, os.path.join(arc_cdir, f"px_{y}.{ext}"))
                elif ext == "png":
                    q = os.path.join(arc_cdir, f"px_{y}.png")
                    if os.path.exists(q):
                        os.remove(q)
            log(f"{cid}: 画素 {y}年（公開中のものを写した）")
            made += 1; continue
        src = get_src()
        dates = src.clear_dates(cell["bbox"], win[0], (we + dt.timedelta(days=1)).isoformat(), MIN_CLEAR)
        layers = {}
        if dates:
            paddy = cached_paddy(src, cdir, fc, g)
            n_paddy = max(1, int(paddy.sum()))
            for d, a in src.ndvi(cell["bbox"], dates, g).items():
                a = a * paddy
                if (a > 0).sum() / n_paddy >= MIN_CLEAR:
                    layers[d] = a
        n = put_year(arc_cdir, y, win, mask, g, dates, layers)
        log(f"{cid}: 画素 {y}年（作った: 晴れた日 {len(dates)}日のうち保存 {n}日）")
        made += 1
    return made


def need_cell(data_dir, out_dir, cell, years, mask, cfg):
    cid = cell["id"]
    parcels = read_json(os.path.join(data_dir, "cells", cid, "parcels.geojson"))
    if not parcels or not parcels.get("features"):
        return False
    g = grid_for(cell["bbox"])
    for y in years:
        ws, we = year_window(cfg, y)
        if not arc_done(os.path.join(out_dir, "cells", cid), y, mask, g, [ws.isoformat(), we.isoformat()]):
            return True
    return False


def main():
    import argparse, sys
    from common import load_config
    from update_ndvi import mask_id
    from trend import run_all
    ap = argparse.ArgumentParser(description="古い年の画素データを px-archive 用に作る")
    ap.add_argument("--data", required=True, help="公開中のデータ（gh-pages の data）")
    ap.add_argument("--out", required=True, help="置き場（px-archive ブランチを取り出した所）")
    ap.add_argument("--years", help="例: 2022-2025 または 2022,2024（省略すると、確定した年すべて）")
    ap.add_argument("--minutes", type=float, default=300, help="使う時間の上限（分）")
    ap.add_argument("--deadline", type=float, help="この時刻（UNIX 秒）までに終える")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--backend", choices=["gee", "fake"], default="gee")
    ap.add_argument("--cells", help="カンマ区切りでセルIDを限定（テスト用）")
    ap.add_argument("--today", help="今日の日付（YYYY-MM-DD。テスト用）")
    args = ap.parse_args()
    def log(m): print(time.strftime("%H:%M:%S"), m, flush=True)
    cfg = load_config()
    index = read_json(os.path.join(args.data, "index.json"))
    if not index:
        sys.exit("index.json がありません")
    today = dt.date.fromisoformat(args.today) if args.today else dt.datetime.now(dt.timezone.utc).date()
    years = parse_years(args.years, cfg, today); mask = mask_id(cfg)
    want = set(args.cells.split(",")) if args.cells else None
    cells = [c for c in index["cells"] if want is None or c["id"] in want]
    limit = args.minutes * 60
    if args.deadline:
        limit = min(limit, args.deadline - time.time())
    lock = threading.Lock(); box = {}
    def get_src():                                  # Earth Engine は、作る必要が出たときに1回だけ始める
        with lock:
            if "src" not in box:
                if args.backend == "gee":
                    from update_ndvi import GEEBackend
                    box["src"] = GEEArc(GEEBackend(cfg))
                else:
                    box["src"] = FakeArc()
            return box["src"]
    need = lambda c: need_cell(args.data, args.out, c, years, mask, cfg)
    log(f"年: {years[0] if years else '-'}〜{years[-1] if years else '-'}（{len(years)}年）, セル {len(cells)}")
    made = 0
    if years and limit > 60:
        made = run_all(None, cfg, args.data, cells, today, mask, log, limit, args.workers,
                       task=lambda c, dl: update_cell(get_src, cfg, args.data, args.out, c, years, mask, today, log, dl),
                       need=need, name="画素の置き場")
    left = sum(1 for c in cells if need(c))
    if years and not os.path.exists(os.path.join(args.out, "README.md")):
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "README.md"), "w", encoding="utf-8") as f:
            f.write("圃場内マップの画素データ（古い年）。pipeline/px_archive.py が作る。画面は raw.githubusercontent.com からここを読む。\n")
    log(f"作った（写した）: {made}セル・年, まだの セル: {left}")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            f.write(f"made={made}\nleft={left}\n")


if __name__ == "__main__":
    main()
