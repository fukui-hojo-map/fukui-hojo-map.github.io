"""高低差マップ（絶対表示）の「周りの田との高さの差」を、比べる範囲（半径）ごとに作る。

入力 : cells/<id>/parcels.geojson（田の形） と cells/<id>/wet.json の terr.p[pid][0]（田の標高の中央値 cm。pipeline/wet.py が国土地理院の標高タイルから作る）
出力 : cells/<id>/elev.json  {"v":1,"r":[1000,2000],"y":["R01","H25",…],"p":{pid:[差cm（半径1000m）, 差cm（半径2000m）, 測量年度の番号]}}   差 = その田の標高 − 周り（半径内の、ほかの田）の標高の中央値。プラスが高い。
       周りの田が8枚に満たないときは null。測量年度の番号は y の何番目か（5mレーザー測量の標高の田だけ。ほかは null）
測量年度: 国土地理院の「5mメッシュDEM（航空レーザ測量）提供範囲と作業年度」(fgd_dem5a_area_dtil, z12 の geojson・場所ごとに原典資料の作業年度 H19〜R06 が違う・重ならない)を田の位置で引く。
       月は公開データにない（年度だけ）。タイルが取れなかったときは、前に作った elev.json の年度を引き継ぐ
半径 300m の差は、wet.json の terr.p[pid][1]（田の縁から300m・10m画素の中央値）をそのまま使う（ここでは作らない）。
全セルの田をまとめて探す（セルの境目で切れない）。田の位置は形の頂点の平均、周りは田の位置どうしの距離（地球を球とした弦の長さ）で数える。
使い方: python pipeline/elev_rel.py --data site/data [--out site/data] [--area-dir 試験用のタイルの置き場]     （数分。Earth Engine は使わない）
"""
import os, sys, json, math, argparse, time, urllib.request, urllib.error
import numpy as np

AREA_URL = "https://cyberjapandata.gsi.go.jp/xyz/fgd_dem5a_area_dtil/12/{x}/{y}.geojson"
AREA_Z = 12
UA = "Mozilla/5.0 (fukui-hojo-map elev_rel; +https://fukui-hojo-map.github.io)"
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
            rows.append((cid, pid, sum(p[0] for p in ring) / len(ring), sum(p[1] for p in ring) / len(ring), a[0], a[8] if len(a) > 8 else None))
    return rows


def area_tile(x, y, area_dir=None):
    """z12 の1枚の作業年度 geojson の特徴 [(作成年度, geometry)]。5mレーザーがない（404）・試験の置き場にない = []。取れなかった = None"""
    if area_dir:
        f = os.path.join(area_dir, f"{x}_{y}.geojson")
        if not os.path.exists(f):
            return []
        return [(ft["properties"].get("作成年度"), ft["geometry"]) for ft in json.load(open(f, encoding="utf-8"))["features"]]
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(AREA_URL.format(x=x, y=y), headers={"User-Agent": UA}), timeout=60) as r:
                return [(ft["properties"].get("作成年度"), ft["geometry"]) for ft in json.loads(r.read().decode("utf-8"))["features"]]
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return []
        except Exception:
            pass
        time.sleep(2 * (k + 1))
    return None


def survey_years(lon, lat, area_dir=None):
    """各点の測量年度（'R01' など）。5mレーザーの範囲の外は None、タイルが取れなかったものは False。z12 のタイルごとに、重ならない範囲の中に入るかを見る"""
    import shapely
    from shapely.geometry import shape
    n = 2 ** AREA_Z; out = np.full(len(lon), None, dtype=object)
    tx = np.floor((lon + 180) / 360 * n).astype(int)
    ty = np.floor((1 - np.log(np.tan(np.radians(lat)) + 1 / np.cos(np.radians(lat))) / math.pi) / 2 * n).astype(int)
    keys = sorted(set(zip(tx.tolist(), ty.tolist())))
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=1 if area_dir else 6) as ex:   # タイルは100枚あまり。6本ずつ取る
        tiles = list(ex.map(lambda k: area_tile(k[0], k[1], area_dir), keys))
    for (x, y), feats in zip(keys, tiles):
        sel = np.where((tx == x) & (ty == y))[0]
        if feats is None:
            out[sel] = False; print(f"作業年度のタイル {x}/{y} が取れません（前の年度を引き継ぎます）", file=sys.stderr); continue
        for code, g in feats:
            hit = shapely.contains_xy(shape(g).buffer(0), lon[sel], lat[sel])
            for j in sel[hit]:
                out[j] = code
    return out


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
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--out"); ap.add_argument("--area-dir"); a = ap.parse_args()
    out_dir = a.out or a.data
    rows = load(a.data)
    if not rows:
        print("標高のある田がありません"); return 1
    lon = np.array([r[2] for r in rows]); lat = np.array([r[3] for r in rows]); elev = np.array([r[4] for r in rows], float)
    rel = rel_all(lon, lat, elev)
    laser = [i for i, r in enumerate(rows) if r[5] == 0]   # 標高が5mレーザー測量の田だけ年度を引く
    yrs = [None] * len(rows)
    if laser:
        got = survey_years(lon[laser], lat[laser], a.area_dir)
        for k, i in enumerate(laser):
            yrs[i] = got[k]
    per = {}; idx = {}
    for i, r in enumerate(rows):
        per.setdefault(r[0], {})[r[1]] = [None if np.isnan(rel[R][i]) else int(round(rel[R][i])) for R in RADII]
        idx.setdefault(r[0], []).append(i)
    n = 0; nyear = 0; nlost = 0
    for cid, p in per.items():
        path = os.path.join(out_dir, "cells", cid, "elev.json"); os.makedirs(os.path.dirname(path), exist_ok=True)
        prev = {}
        if os.path.exists(path):
            try:
                o = json.load(open(path, encoding="utf-8")); prev = {pid: (o["y"][q[2]] if len(q) > 2 and q[2] is not None else None) for pid, q in o["p"].items()} if o.get("y") else {}
            except Exception:
                prev = {}
        codes = []
        for i in idx[cid]:
            r = rows[i]; c = yrs[i]
            if c is False:
                c = prev.get(r[1]); nlost += 1
            if c is not None:
                if c not in codes: codes.append(c)
                per[cid][r[1]].append(codes.index(c)); nyear += 1
            else:
                per[cid][r[1]].append(None)
        body = json.dumps({"v": 1, "r": list(RADII), "y": codes, "p": p}, ensure_ascii=False, separators=(",", ":"))
        if os.path.exists(path) and open(path, encoding="utf-8").read() == body:
            continue
        open(path, "w", encoding="utf-8").write(body); n += 1
    print(f"田 {len(rows)} 枚・セル {len(per)} 枚（書き直し {n}）。半径 {list(RADII)}m。測量年度が付いた田 {nyear} 枚（レーザー測量の標高の田 {len(laser)} 枚・タイルが取れず前の年度で補ったもの {nlost}）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
