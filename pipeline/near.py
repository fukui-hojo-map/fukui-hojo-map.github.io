"""住宅地の近さと田の中の鉄塔（田ごと）を作る。圃場スコアの「住宅地の近さ」の減点と、作業効率性（鉄塔の回り込み）に使う。

cells/<id>/near.json: {"v":1, "src":"worldcover-v200+osm", "tw": 鉄塔のデータがあれば 1, "p":{"<pid>":[b30, b100, towers]}}
  b30 / b100 = 田の縁から 0〜30m / 30〜100m の輪の中の、市街地・建物（ESA WorldCover 2021 の Built-up = 50）の割合（%）
  towers     = 田の中（縁から内側3mまで含む）に立つ送電線の鉄塔の数（OpenStreetMap の power=tower）
鉄塔は県の範囲を Overpass API から1回だけ取り、data/osm_towers.json に置いて使い回す（取れなかった年は tw=0 のまま、次の回にもう一度取る）
"""
import os, math, json, urllib.request, urllib.parse
from common import write_json, read_json, with_retry
from pixels import grid_for

VERSION = 1
PAD_DEG = (0.0015, 0.0012)       # セルの外に広げる幅（経度, 緯度）: 約130m
RINGS_M = (30, 100)
FUKUI = (35.30, 135.40, 36.35, 136.90)   # 南, 西, 北, 東
OVERPASS = ("https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter")


class GEENear:
    def __init__(self, backend, px):
        self.ee = backend.ee; self.px = px

    def built(self, g):
        import numpy as np
        ee = self.ee
        img = ee.ImageCollection("ESA/WorldCover/v200").first().select("Map").eq(50).unmask(0).byte().rename("b")
        return np.asarray(with_retry(lambda: ee.data.computePixels(self.px._req(img, g)))["b"], dtype=np.uint8)


class FakeNear:
    """動作確認用: セルの東の端から100mは市街地"""
    def built(self, g):
        import numpy as np
        a = np.zeros((g["h"], g["w"]), np.uint8); a[:, -int(100 / 10) - int(PAD_DEG[0] * 111000 * 0.81 / 10):] = 1
        return a


def towers(data_dir, log=print):
    """県の鉄塔 [[lon, lat], ...]。取れなければ None"""
    path = os.path.join(data_dir, "osm_towers.json")
    t = read_json(path)
    if t and t.get("pts") is not None:
        return t["pts"]
    s, w, n, e = FUKUI
    q = f'[out:json][timeout:120];node["power"="tower"]({s},{w},{n},{e});out;'
    for url in OVERPASS:
        try:
            req = urllib.request.Request(url, data=urllib.parse.urlencode({"data": q}).encode(), headers={"User-Agent": "tanbo-ndvi"})
            with urllib.request.urlopen(req, timeout=180) as r:
                js = json.load(r)
            pts = [[round(el["lon"], 6), round(el["lat"], 6)] for el in js.get("elements", []) if el.get("type") == "node"]
            write_json(path, {"src": "OpenStreetMap power=tower", "pts": pts})
            log(f"鉄塔 {len(pts)} 基（OpenStreetMap）")
            return pts
        except Exception as ex:
            log(f"鉄塔の取得に失敗 {url}: {ex}")
    return None


def update_cell(src, cfg, data_dir, index, cell, log, deadline=None, tw=None):
    """1セルの near.json を作る。作ったら 1 を返す"""
    import numpy as np, shapely
    from shapely.geometry import shape
    cdir = os.path.join(data_dir, "cells", cell["id"])
    parcels = read_json(os.path.join(cdir, "parcels.geojson"))
    if not parcels or not parcels.get("features"):
        return 0
    w, s, e, n = cell["bbox"]
    g = grid_for([w - PAD_DEG[0], s - PAD_DEG[1], e + PAD_DEG[0], n + PAD_DEG[1]])   # 格子1つ ≒ 地上10m
    def to_grid(xy):
        x = 6378137.0 * np.radians(xy[:, 0]); y = 6378137.0 * np.log(np.tan(np.pi / 4 + np.radians(xy[:, 1]) / 2))
        return np.column_stack([(x - g["x0"]) / g["px"], (g["y1"] - y) / g["px"]])
    built = src.built(g).astype(bool)
    yy, xx = np.mgrid[0:g["h"], 0:g["w"]]; cx = xx + 0.5; cy = yy + 0.5
    tp = None
    if tw is not None:
        a = np.array([p for p in tw if w - 0.001 <= p[0] <= e + 0.001 and s - 0.001 <= p[1] <= n + 0.001], float).reshape(-1, 2)
        tp = shapely.points(to_grid(a)) if len(a) else []
    out = {}
    for f in parcels["features"]:
        pid = f["properties"]["pid"]
        try:
            gg = shapely.transform(shape(f["geometry"]), to_grid).buffer(0)
            far = gg.buffer(RINGS_M[-1] / 10.0)
            x0, y0, x1, y1 = far.bounds
            c0, c1, r0, r1 = max(0, math.floor(x0)), min(g["w"], math.ceil(x1)), max(0, math.floor(y0)), min(g["h"], math.ceil(y1))
            X, Y = cx[r0:r1, c0:c1], cy[r0:r1, c0:c1]; B = built[r0:r1, c0:c1]
            inside = shapely.contains_xy(gg, X, Y)
            dist = np.asarray(shapely.distance(gg, shapely.points(X.ravel(), Y.ravel()))).reshape(X.shape) * 10.0
            dist[inside] = 0
            row = []; lo = 0
            for hi in RINGS_M:
                ring = (dist > lo) & (dist <= hi) & ~inside
                row.append(round(100 * B[ring].mean()) if ring.any() else 0); lo = hi
            nt = 0
            if tp is not None and len(tp):
                gf = shapely.Polygon(gg.exterior) if gg.geom_type == "Polygon" else shapely.MultiPolygon([shapely.Polygon(q.exterior) for q in gg.geoms])   # 筆ポリゴンは鉄塔の足元を穴で抜いてあるので、穴を埋めて数える
                gi = gf.buffer(-0.3)
                nt = int(shapely.contains(gi if gi.area > 0 else gf, tp).sum())   # 縁から内側 3m より中（畦の上の鉄塔は数えない）
            row.append(nt)
            out[pid] = row
        except Exception:
            continue
    write_json(os.path.join(cdir, "near.json"), {"v": VERSION, "src": "worldcover-v200+osm", "tw": 1 if tw is not None else 0, "p": out})
    return 1


def need(data_dir, cell, tw_ok=True):
    m = read_json(os.path.join(data_dir, "cells", cell["id"], "near.json")) or {}
    return m.get("v") != VERSION or (tw_ok and not m.get("tw"))
