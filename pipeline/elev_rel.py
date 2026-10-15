"""高低差マップ（絶対表示）の「周りの田との高さの差」を、比べる範囲（半径）ごとに作る。

入力 : cells/<id>/parcels.geojson（田の形） と cells/<id>/wet.json の terr.p[pid][0]（田の標高の中央値 cm。pipeline/wet.py が国土地理院の標高タイルから作る）
出力 : cells/<id>/elev.json  {"v":1,"r":[1000,2000],"p":{pid:[差cm（半径1000m）, 差cm（半径2000m）]}}   差 = その田の標高 − 周り（半径内の、ほかの田）の標高の中央値。プラスが高い。
       周りの田が8枚に満たないときは null
半径 300m の差は、wet.json の terr.p[pid][1]（田の縁から300m・10m画素の中央値）をそのまま使う（ここでは作らない）。
全セルの田をまとめて探す（セルの境目で切れない）。田の位置は形の頂点の平均、周りは田の位置どうしの距離（地球を球とした弦の長さ）で数える。
使い方: python pipeline/elev_rel.py --data site/data [--out site/data]     （1〜2分。Earth Engine は使わない）
"""
import os, sys, json, math, argparse
import numpy as np

RADII = (1000, 2000)
MIN_NEIGHBORS = 8
CHUNK = 3000


def load(data):
    idx = json.load(open(os.path.join(data, "index.json"), encoding="utf-8"))
    rows = []   # (cell, pid, lon, lat, 標高cm)
    for c in idx["cells"]:
        cid = c["id"]; cdir = os.path.join(data, "cells", cid)
        pf, wf = os.path.join(cdir, "parcels.geojson"), os.path.join(cdir, "wet.json")
        if not (os.path.exists(pf) and os.path.exists(wf)):
            continue
        try:
            t = ((json.load(open(wf, encoding="utf-8")).get("terr") or {}).get("p")) or {}
            feats = json.load(open(pf, encoding="utf-8"))["features"]
        except Exception as e:
            print(f"{cid}: 読めません（{e}）", file=sys.stderr); continue
        for f in feats:
            pid = f["properties"].get("pid"); a = t.get(pid)
            if not a or a[0] is None:
                continue
            g = f["geometry"]; ring = g["coordinates"][0] if g["type"] == "Polygon" else g["coordinates"][0][0]
            rows.append((cid, pid, sum(p[0] for p in ring) / len(ring), sum(p[1] for p in ring) / len(ring), a[0]))
    return rows


def rel_all(lon, lat, elev, radii=RADII):
    from scipy.spatial import cKDTree
    la, lo = np.radians(lat), np.radians(lon); Rm = 6371000.0   # 地球の中心を原点にした3次元の座標（弦の長さ ≒ 地表の距離）。緯度による経度1度の長さの違いが出ない
    xy = np.column_stack([Rm * np.cos(la) * np.cos(lo), Rm * np.cos(la) * np.sin(lo), Rm * np.sin(la)]); tree = cKDTree(xy)
    out = {}
    for R in radii:
        rel = np.full(len(elev), np.nan)
        for s in range(0, len(elev), CHUNK):
            for k, nb in enumerate(tree.query_ball_point(xy[s:s + CHUNK], R)):
                i = s + k; nb = [j for j in nb if j != i]
                if len(nb) >= MIN_NEIGHBORS:
                    rel[i] = elev[i] - np.median(elev[nb])
        out[R] = rel
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--out"); a = ap.parse_args()
    out_dir = a.out or a.data
    rows = load(a.data)
    if not rows:
        print("標高のある田がありません"); return 1
    lon = np.array([r[2] for r in rows]); lat = np.array([r[3] for r in rows]); elev = np.array([r[4] for r in rows], float)
    rel = rel_all(lon, lat, elev)
    per = {}
    for i, r in enumerate(rows):
        per.setdefault(r[0], {})[r[1]] = [None if np.isnan(rel[R][i]) else int(round(rel[R][i])) for R in RADII]
    n = 0
    for cid, p in per.items():
        path = os.path.join(out_dir, "cells", cid, "elev.json"); os.makedirs(os.path.dirname(path), exist_ok=True)
        body = json.dumps({"v": 1, "r": list(RADII), "p": p}, ensure_ascii=False, separators=(",", ":"))
        if os.path.exists(path) and open(path, encoding="utf-8").read() == body:
            continue
        open(path, "w", encoding="utf-8").write(body); n += 1
    print(f"田 {len(rows)} 枚・セル {len(per)} 枚（書き直し {n}）。半径 {list(RADII)}m")
    return 0


if __name__ == "__main__":
    sys.exit(main())
