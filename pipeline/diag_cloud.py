"""T19 診断（探索）: 積雲（点々とある雲）の近くの圃場で、Sentinel-2 の NDVI が下がるか。
圃場（inner.geojson の内側の形）ごと・日ごとに、今の判定の NDVI と各バンド、Cloud Score+、AOT、
雲（cs_cdf<0.6 か SCL 8・9・10）と雲の影（SCL 3）までの距離、近く（1km・3km）の雲の割合を出す。
同じ日の画像を別の方法で補正した NASA HLS S30（LaSRC・Fmask。雲のそばの画素の印がある）と、Landsat（HLS L30）、
MODIS（Terra・Aqua、毎日、250m）の、稲の田の上の NDVI も出す。地図の画像（真色・NDVI・雲までの距離）も作る。値の補正はしない。
使い方: python pipeline/diag_cloud.py --data <公開データの data> --out <出力> --cells 5449-5452x2161-2164 [--dates 2026-07-08,2026-07-23]
"""
import os, sys, json, glob, time, argparse, gzip, csv, datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed

CHUNK = 1500

def init_ee():
    import ee
    sa, key, project = os.environ.get("EE_SERVICE_ACCOUNT"), os.environ.get("EE_PRIVATE_KEY"), os.environ.get("EE_PROJECT")
    if sa and key: ee.Initialize(ee.ServiceAccountCredentials(sa, key_data=key), project=project)
    else: ee.Initialize(project=project)
    return ee

def parse_cells(spec):
    out = []
    for part in spec.split(","):
        xs, ys = part.split("x"); x0, x1 = (int(v) for v in xs.split("-")); y0, y1 = (int(v) for v in ys.split("-"))
        out += [f"{x}_{y}" for x in range(x0, x1 + 1) for y in range(y0, y1 + 1)]
    return out

def load_fields(data, cells, keep, only=None):
    """稲（2025 か 2026 に確信度80以上）の田と、keep の田。[(pid, geometry, rice25, rice26)]。only（セル→田の集合）があれば、その田だけ"""
    feats, bbox = [], [999, 999, -999, -999]
    for c in cells:
        fi, fc = os.path.join(data, "cells", c, "inner.geojson"), os.path.join(data, "cells", c, "crop_s.json")
        if not (os.path.exists(fi) and os.path.exists(fc)): continue
        cs = json.load(open(fc, encoding="utf-8")); yrs = cs["years"]
        for f in json.load(open(fi, encoding="utf-8"))["features"]:
            pid = f["properties"]["pid"]; v = cs["p"].get(pid)
            r = lambda y: bool(v) and y in yrs and v[yrs.index(y)][0] == "稲" and v[yrs.index(y)][1] >= 80
            r25, r26 = r(2025), r(2026)
            if only is not None:
                if pid not in only.get(c, ()): continue
            elif not (r25 or r26 or pid in keep): continue
            feats.append((pid, f["geometry"], r25, r26))
            for x, y in _coords(f["geometry"]):
                bbox = [min(bbox[0], x), min(bbox[1], y), max(bbox[2], x), max(bbox[3], y)]
    return feats, bbox

def _coords(g):
    def walk(a):
        if isinstance(a[0], (int, float)): yield a
        else:
            for b in a: yield from walk(b)
    yield from walk(g["coordinates"])

def retry(fn, tries=4, what=""):
    last = None
    for k in range(tries):
        try: return fn()
        except Exception as e:
            last = e; print("やり直し", what, str(e)[:200], flush=True); time.sleep(8 * (k + 1))
    raise last

# ---------- Sentinel-2（今の判定と同じ: S2_SR_HARMONIZED + Cloud Score+） ----------
def s2_dates(ee, region, start, end):
    sr = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(region).filterDate(start, end)
    cp = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED").filterBounds(region).filterDate(start, end)
    sr = sr.filter(ee.Filter.inList("system:index", cp.aggregate_array("system:index")))
    info = retry(lambda: sr.map(lambda i: ee.Feature(None, {"t": i.get("system:time_start"), "sc": i.get("SPACECRAFT_NAME"), "o": i.get("SENSING_ORBIT_NUMBER"),
        "cl": i.get("CLOUDY_PIXEL_PERCENTAGE"), "sza": i.get("MEAN_SOLAR_ZENITH_ANGLE"), "vza": i.get("MEAN_INCIDENCE_ZENITH_ANGLE_B8"), "tile": i.get("MGRS_TILE"),
        "aotm": i.get("AOT_RETRIEVAL_METHOD"), "pb": i.get("PROCESSING_BASELINE")})).getInfo(), what="s2 一覧")
    meta = {}
    for f in info["features"]:
        p = f["properties"]; d = dt.datetime.fromtimestamp(p["t"] / 1000, dt.timezone.utc).strftime("%Y-%m-%d")
        meta.setdefault(d, []).append({k: p.get(k) for k in ("sc", "o", "cl", "sza", "vza", "tile", "aotm", "pb")})
    return meta

def s2_stack(ee, d, region):
    proj = ee.Projection("EPSG:32653").atScale(10)
    d0 = ee.Date(d); d1 = d0.advance(1, "day")
    sr = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(region).filterDate(d0, d1)
    cp = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED").filterBounds(region).filterDate(d0, d1)
    sr = sr.filter(ee.Filter.inList("system:index", cp.aggregate_array("system:index")))
    img = sr.linkCollection(cp, ["cs", "cs_cdf"]).mosaic().setDefaultProjection(proj)
    scl, cdf, cs = img.select("SCL"), img.select("cs_cdf"), img.select("cs")
    ok0 = cdf.gte(0.6).And(scl.neq(0)).And(scl.neq(1)).And(scl.neq(11))
    ndvi = img.normalizedDifference(["B8", "B4"]).rename("ndvi")
    refl = img.select(["B2", "B3", "B4", "B8", "B11"]).multiply(1e-4).rename(["b2", "b3", "b4", "b8", "b11"])
    cloud = cdf.lt(0.6).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10))
    cloud_scl = scl.eq(8).Or(scl.eq(9)).Or(scl.eq(10))
    shadow = scl.eq(3)
    def frac(m, s): return m.unmask(0).reproject(proj).reduceResolution(ee.Reducer.mean(), False, 1024).reproject(proj.atScale(s))
    def dist(m):    # 30m で 25%以上が雲の画素までの距離（m、最大4.5km）
        src = frac(m, 30).gte(0.25)
        return src.fastDistanceTransform(150, "pixels", "squared_euclidean").sqrt().multiply(30).reproject(proj.atScale(30)).unmask(4500).min(4500)
    f1k = frac(cloud, 30).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(1000, "meters")).reproject(proj.atScale(30)).rename("f1k")
    f3k = frac(cloud, 90).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(3000, "meters")).reproject(proj.atScale(90)).rename("f3k")
    one = img.select("B4").gte(-1e9).rename("one")
    b10 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").filterBounds(region).filterDate(d0, d1).mosaic().select("B10").multiply(1e-4).rename("b10")
    stack = (one.addBands(ndvi.updateMask(ok0).rename("ndvi0")).addBands(refl.updateMask(ok0)).addBands(b10)
             .addBands(cs.rename("cs")).addBands(cdf.rename("cdf")).addBands(img.select("AOT").multiply(0.001).rename("aot"))
             .addBands(dist(cloud).rename("dc")).addBands(dist(cloud_scl).rename("dcs")).addBands(dist(shadow).rename("ds"))
             .addBands(f1k).addBands(f3k))
    return stack, img, ndvi.updateMask(ok0), dist(cloud)

S2_COLS = ["pid", "n_all", "n0", "ndvi0", "b2", "b3", "b4", "b8", "b11", "cs", "cdf", "aot", "dc", "dcs", "ds", "f1k", "f3k", "b10"]

def s2_rows(ee, d, region, fc):
    stack, *_ = s2_stack(ee, d, region)
    red = ee.Reducer.mean().combine(ee.Reducer.count(), "", True)
    res = retry(lambda: stack.reduceRegions(collection=fc, reducer=red, scale=10, tileScale=4).getInfo(), what=f"s2 {d}")
    rows = []
    for f in res["features"]:
        p = f["properties"]
        if not p.get("one_count"): continue
        g = lambda k: p.get(k + "_mean")
        rows.append([p["pid"], p.get("one_count"), p.get("ndvi0_count") or 0] + [None if g(k) is None else round(g(k), 5) for k in S2_COLS[3:]])
    return rows

# ---------- HLS（S30 = 同じ Sentinel-2 の画像を LaSRC で補正、L30 = Landsat） ----------
def hls_dates(ee, coll, region, start, end):
    ts = retry(lambda: ee.ImageCollection(coll).filterBounds(region).filterDate(start, end).aggregate_array("system:time_start").getInfo(), what="hls 一覧")
    return sorted({dt.datetime.fromtimestamp(t / 1000, dt.timezone.utc).strftime("%Y-%m-%d") for t in ts})

def hls_rows(ee, coll, nir, d, region, fc):
    d0 = ee.Date(d); d1 = d0.advance(1, "day")
    img = ee.ImageCollection(coll).filterBounds(region).filterDate(d0, d1).mosaic().setDefaultProjection(ee.Projection("EPSG:32653").atScale(30))
    fm = img.select("Fmask"); bit = lambda b: fm.rightShift(b).bitwiseAnd(1).eq(1); aer = fm.rightShift(6).bitwiseAnd(3)
    ndvi = img.normalizedDifference([nir, "B4"])
    strict = bit(1).Or(bit(2)).Or(bit(3)).Or(bit(4)).Or(aer.eq(3)).Not()
    noadj = bit(1).Or(bit(3)).Or(bit(4)).Or(aer.eq(3)).Not()
    st = (img.select("B4").gte(-1e9).rename("one").addBands(ndvi.updateMask(strict).rename("ns")).addBands(ndvi.updateMask(noadj).rename("na"))
          .addBands(img.select(["B2", "B4", nir]).updateMask(strict).rename(["b2", "b4", "nir"])).addBands(aer.rename("aer")))
    red = ee.Reducer.mean().combine(ee.Reducer.count(), "", True)
    res = retry(lambda: st.reduceRegions(collection=fc, reducer=red, scale=30, tileScale=4).getInfo(), what=f"hls {coll} {d}")
    rows = []
    for f in res["features"]:
        p = f["properties"]
        if not p.get("one_count"): continue
        g = lambda k: None if p.get(k + "_mean") is None else round(p[k + "_mean"], 5)
        rows.append([p["pid"], p.get("one_count"), p.get("ns_count") or 0, g("ns"), p.get("na_count") or 0, g("na"), g("b2"), g("b4"), g("nir"), g("aer")])
    return rows

HLS_COLS = ["pid", "n_all", "n_strict", "ndvi_strict", "n_noadj", "ndvi_noadj", "b2", "b4", "nir", "aer"]

# ---------- MODIS（毎日、250m。稲の田が多い画素だけ） ----------
def modis_rows(ee, rice_fc, region, start, end):
    proj10 = ee.Projection("EPSG:32653").atScale(10)
    gq0 = ee.ImageCollection("MODIS/061/MOD09GQ").filterDate(start, end).first()
    mproj = gq0.select("sur_refl_b01").projection()
    frac = ee.Image(0).byte().paint(rice_fc, 1).reproject(proj10).reduceResolution(ee.Reducer.mean(), False, 1024).reproject(mproj)
    out = []
    days = (dt.date.fromisoformat(end) - dt.date.fromisoformat(start)).days
    for sat in ("MOD", "MYD"):
        gq = ee.ImageCollection(f"MODIS/061/{sat}09GQ").filterDate(start, end)
        ga = ee.ImageCollection(f"MODIS/061/{sat}09GA").filterDate(start, end)
        def per(i):
            i = ee.Image(i); dd = i.date()
            a = ga.filterDate(dd, dd.advance(1, "day")).first()
            stt = a.select("state_1km"); b = lambda s, n=1: stt.rightShift(s).bitwiseAnd((1 << n) - 1)
            clear = b(0, 2).eq(0).And(b(2).eq(0)).And(b(8, 2).eq(0)).And(b(10).eq(0)).And(b(13).eq(0)).And(b(6, 2).neq(3))
            good = i.select("QC_250m").bitwiseAnd(3).eq(0)
            r, n = i.select("sur_refl_b01").multiply(1e-4), i.select("sur_refl_b02").multiply(1e-4)
            nd = n.subtract(r).divide(n.add(r)).rename("nd").updateMask(clear).updateMask(good).updateMask(frac.gte(0.5)).reproject(mproj)
            vz = a.select("SensorZenith").multiply(0.01).rename("vz").reproject(mproj)
            s = nd.addBands(vz.updateMask(nd.mask())).reduceRegion(ee.Reducer.percentile([25, 50, 75]).combine(ee.Reducer.count(), "", True), region, 250, maxPixels=1e8)
            tot = frac.gte(0.5).rename("t").reduceRegion(ee.Reducer.sum(), region, 250, maxPixels=1e8)
            return ee.Feature(None, s).set({"date": dd.format("YYYY-MM-dd"), "sat": sat, "tot": tot.get("t")})
        res = retry(lambda: ee.FeatureCollection(gq.map(per)).getInfo(), what=f"modis {sat}")
        for f in res["features"]:
            p = f["properties"]
            out.append([p.get("date"), sat, p.get("tot"), p.get("nd_count"), p.get("nd_p25"), p.get("nd_p50"), p.get("nd_p75"), p.get("vz_p50")])
    return out

MODIS_COLS = ["date", "sat", "rice_px", "n_clear", "ndvi_p25", "ndvi_p50", "ndvi_p75", "vza_p50"]

def write_gz(path, cols, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols); w.writerows(rows)

def thumbs(ee, out, region_box, s2days, l30days):
    import requests
    reg = ee.Geometry.Rectangle(region_box)
    def save(name, img):
        for k in range(3):
            try:
                url = img.getThumbURL({"region": reg, "dimensions": "720", "format": "png"})
                r = requests.get(url, timeout=180); r.raise_for_status(); open(os.path.join(out, name), "wb").write(r.content); print("画像", name, len(r.content), flush=True); return
            except Exception as e:
                print("画像失敗", name, str(e)[:160], flush=True); time.sleep(5 * (k + 1))
    pal = ["a52a2a", "f5e663", "7fc41c", "0b6b1e"]
    for d in s2days:
        stack, img, nd0, dc = s2_stack(ee, d, reg)
        save(f"img_{d}_s2_true.png", img.select(["B4", "B3", "B2"]).visualize(min=0, max=2500, gamma=1.2))
        save(f"img_{d}_s2_ndvi.png", nd0.visualize(min=0.5, max=0.95, palette=pal))
        save(f"img_{d}_s2_dist.png", dc.visualize(min=0, max=3000, palette=["ff0000", "ffaa00", "ffff66", "ffffff"]))
    for d in l30days:
        d0 = ee.Date(d)
        im = ee.ImageCollection("NASA/HLS/HLSL30/v002").filterBounds(reg).filterDate(d0, d0.advance(1, "day")).mosaic()
        save(f"img_{d}_l30_true.png", im.select(["B4", "B3", "B2"]).visualize(min=0, max=0.25, gamma=1.2))
        save(f"img_{d}_l30_true2500.png", im.select(["B4", "B3", "B2"]).visualize(min=0, max=2500, gamma=1.2))   # 反射率の桁が違うときのため
        fm = im.select("Fmask"); okm = fm.bitwiseAnd(0b11110).eq(0)
        save(f"img_{d}_l30_ndvi.png", im.normalizedDifference(["B5", "B4"]).updateMask(okm).visualize(min=0.5, max=0.95, palette=pal))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--out", required=True); ap.add_argument("--cells", default="")
    ap.add_argument("--dates", default=""); ap.add_argument("--periods", default="2025-06-01:2025-10-01,2026-06-01:2026-10-01")
    ap.add_argument("--keep", default="50f6ec09-1cb9-4d11-8bbf-4c40bb9a759f"); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--deadline", type=float, default=0); ap.add_argument("--parts", default="s2,s30,l30,modis,img")
    ap.add_argument("--img-box", default="136.235,36.035,136.32,36.09"); ap.add_argument("--img-s2", default="2026-07-08,2026-07-18,2026-07-23,2026-08-04,2026-08-07")
    ap.add_argument("--img-l30", default="2026-07-08,2026-07-24")
    ap.add_argument("--pids-file", default="", help="「セルID,田のID」の行の一覧。あれば、その田だけ（--cells は不要）")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ee = init_ee(); parts = set(a.parts.split(","))
    only = None
    if a.pids_file:
        only = {}
        for ln in open(a.pids_file, encoding="utf-8"):
            ln = ln.strip()
            if ln and "," in ln:
                c, p = ln.split(",", 1); only.setdefault(c, set()).add(p)
        cells = sorted(only)
    else:
        cells = parse_cells(a.cells)
    feats, bbox = load_fields(a.data, cells, set(a.keep.split(",")), only)
    print(f"セル {len(cells)}、田 {len(feats)}、範囲 {bbox}", flush=True)
    cur = {"cells": cells, "fields": [[p, r25, r26] for p, _, r25, r26 in feats], "bbox": bbox}
    fpath = os.path.join(a.out, "fields.json"); old = None
    if os.path.exists(fpath):
        try: old = json.load(open(fpath))
        except Exception: old = None
    if old is not None and old.get("fields") != cur["fields"]:     # 前回と田の組が違う: 前回の結果は使えない
        for f in glob.glob(os.path.join(a.out, "*.csv.gz")) + glob.glob(os.path.join(a.out, "img_*.png")) + glob.glob(os.path.join(a.out, "s2_meta.json")):
            os.remove(f)
        print("前回の結果は、田の組が違うので消した", flush=True)
    json.dump(cur, open(fpath, "w"))
    region = ee.Geometry.Rectangle(bbox)
    chunks = [ee.FeatureCollection([ee.Feature(g, {"pid": p}) for p, g, _, _ in feats[i:i + CHUNK]]) for i in range(0, len(feats), CHUNK)]
    periods = [p.split(":") for p in a.periods.split(",")]
    late = lambda: a.deadline and time.time() > a.deadline
    jobs = []
    if "s2" in parts:
        meta = {}
        for s, e in periods: meta.update(s2_dates(ee, region, s, e))
        json.dump(meta, open(os.path.join(a.out, "s2_meta.json"), "w"), ensure_ascii=False)
        ds = [d for d in sorted(meta) if not a.dates or d in a.dates.split(",")]
        print("S2 の日", len(ds), flush=True)
        jobs += [("s2", d) for d in ds]
    for src, coll, nir in (("s30", "NASA/HLS/HLSS30/v002", "B8A"), ("l30", "NASA/HLS/HLSL30/v002", "B5")):
        if src not in parts: continue
        ds = []
        for s, e in periods: ds += hls_dates(ee, coll, region, s, e)
        ds = [d for d in ds if not a.dates or d in a.dates.split(",") or src == "l30"]
        print(src, "の日", len(ds), flush=True)
        jobs += [(src, d) for d in ds]
    def run(job):
        src, d = job; path = os.path.join(a.out, f"{src}_{d}.csv.gz")
        if os.path.exists(path): return job, "あり"
        if late(): return job, "時間切れ"
        rows = []
        for fc in chunks:
            rows += s2_rows(ee, d, region, fc) if src == "s2" else hls_rows(ee, {"s30": "NASA/HLS/HLSS30/v002", "l30": "NASA/HLS/HLSL30/v002"}[src], {"s30": "B8A", "l30": "B5"}[src], d, region, fc)
        write_gz(path, S2_COLS if src == "s2" else HLS_COLS, rows)
        return job, len(rows)
    t0 = time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        for fu in as_completed([ex.submit(run, j) for j in jobs]):
            try: job, n = fu.result(); print("済", job, n, f"{time.time() - t0:.0f}s", flush=True)
            except Exception as e: print("失敗", str(e)[:300], flush=True)
    if "modis" in parts and not late():
        rice_fc = ee.FeatureCollection([ee.Feature(g, {}) for p, g, r25, r26 in feats if r25 and r26])
        rows = []
        for s, e in periods:
            s2_ = max(s, s[:4] + "-06-15"); e2 = min(e, e[:4] + "-09-01")
            try: rows += modis_rows(ee, rice_fc, region, s2_, e2)
            except Exception as ex_: print("MODIS 失敗", str(ex_)[:300], flush=True)
        write_gz(os.path.join(a.out, "modis.csv.gz"), MODIS_COLS, rows); print("MODIS", len(rows), flush=True)
    if "img" in parts and not late():
        try: thumbs(ee, a.out, [float(v) for v in a.img_box.split(",")], a.img_s2.split(","), a.img_l30.split(","))
        except Exception as e: print("画像 失敗", str(e)[:300], flush=True)
    print("完了", flush=True)

if __name__ == "__main__": main()
