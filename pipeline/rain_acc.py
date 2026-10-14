# 圃場ごとの雨量（12・24・48時間）を、気象庁の解析雨量から 1km 四方ごとに組み立てる（Actions で 30分ごと）
#
# 気象庁の「今後の雨」のデータ（www.jma.go.jp/bosai/jmatile/data/rasrf/）から、次の2つを使う:
#   ・1km 四方の色の画像（z10 のタイル）: 1時間・3時間・24時間の雨量の階級（1時間は 1mm 未満/1/5/10/20/30/50/80mm）
#   ・地点の値（data.geojson）: 5km 四方（1km 格子 5列×6行）の中の「いちばん多い 1km 格子」の値（1時間は 1mm 刻み、1mm 未満は 0.4）
# 1時間ごとに、1km 格子の値を「階級の下から2割」で置き、まわりとなめらかにつなぎ（階級の幅と 5km の最大の中で）、
# さらに 3時間・24時間の階級と 5km の最大に合うように直す。それを 12・24・48時間で足す。
# 結果は acc.json（1km 格子ごとの 0.1mm 単位）。アメダスと突き合わせた答え合わせ（check.jsonl）も残す。
# 気象庁ホームページのデータを加工して作成（出典の表示は画面で行う）。気象庁の解析（観測）だけで、予報は作らない。
import argparse, base64, io, json, math, os, re, sys, time, threading, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone, timedelta
import numpy as np
from PIL import Image

UA = "Mozilla/5.0 (fukui-hojo-map rain; +https://fukui-hojo-map.github.io)"
J = "https://www.jma.go.jp"
D = J + "/bosai/jmatile/data/rasrf/"
A = J + "/bosai/amedas/"
# 1km 格子（解析雨量の格子: 経度 118° から 45秒＝0.0125° ごと、緯度 48° から 30秒＝1/120° ごと）。福井県とまわり
P0, P1 = 1388, 1516        # 経度 135.35〜136.95
Q0, Q1 = 1386, 1524        # 緯度 36.45〜35.30（北→南）
NX, NY = P1 - P0, Q1 - Q0
TX, TY = range(897, 902), range(400, 405)   # z10 のタイル
CI0, CI1, CJ0, CJ1 = P0 // 5, (P1 - 1) // 5, Q0 // 6, (Q1 - 1) // 6   # 5km 四方の番号の範囲
COL = [(242, 242, 255), (160, 210, 255), (33, 140, 255), (0, 65, 255), (250, 245, 0), (255, 153, 0), (255, 40, 0), (180, 0, 104)]
# 階級の範囲（mm）。1時間は 1mm 刻みの整数（1mm 未満は別）。3・24時間は 0.2mm 刻み
B1 = [(0, 0), (1, 4), (5, 9), (10, 19), (20, 29), (30, 49), (50, 79), (80, None)]
B3 = [(1, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 150), (150, None)]
B24 = [(1, 50), (50, 80), (80, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, None)]
FRAC = 0.2       # 階級の下から2割（雨量計で確かめた。2〜2.5割で、ずれがいちばん小さい）
C0 = 0.12        # 1mm 未満の1時間（雨量計の平均）。気象庁の積算では 0.4 と数える
SMOOTH = 2
WINS = (12, 24, 48)
KEEP_H = 74      # 1時間ごとの記録をとっておく長さ（48時間＋24時間の制約）

def log(*a): print(*a, flush=True)
_lk = threading.Lock()
def get(url, timeout=25):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for k in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                b = r.read()
                if b[:2] == b"\x1f\x8b":
                    import gzip; b = gzip.decompress(b)
                return r.status, b
        except urllib.error.HTTPError as e:
            if e.code in (404, 403, 410): return e.code, b""
            if k == 2: return e.code, b""
        except Exception as e:
            if k == 2: return None, str(e).encode()
        time.sleep(1.5 * (k + 1))

# ---- 1km 格子の中心が、どのタイルのどの画素か ----
def _mesh_pix():
    n = 1024; tix = np.zeros((NY, NX), int); pxs = np.zeros((NY, NX), int); pys = np.zeros((NY, NX), int)
    tl = [(x, y) for x in TX for y in TY]
    for qi in range(NY):
        lat = 48 - (Q0 + qi + 0.5) / 120; yf = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n; y = int(yf)
        for pi in range(NX):
            lon = 118 + (P0 + pi + 0.5) * 0.0125; xf = (lon + 180) / 360 * n; x = int(xf)
            tix[qi, pi] = tl.index((x, y)); pxs[qi, pi] = int((xf - x) * 256); pys[qi, pi] = int((yf - y) * 256)
    return tl, tix, pxs, pys
TL, TIX, PXS, PYS = _mesh_pix()
CELLI = (np.arange(P0, P1) // 5 - CI0)[None, :].repeat(NY, 0)
CELLJ = (np.arange(Q0, Q1) // 6 - CJ0)[:, None].repeat(NX, 1)

ND_PTS = [(lat, lon) for lat in (35.45, 35.65, 35.85, 36.05, 36.25) for lon in (135.6, 135.9, 136.2, 136.5, 136.75)]   # 福井の陸のあたりの点
def _in_rings(lat, lon, rings):
    """多角形（外周と穴）の中か（偶奇の規則）"""
    inside = False
    for ring in rings:
        n = len(ring)
        for i in range(n):
            x1, y1 = ring[i][0], ring[i][1]; x2, y2 = ring[(i + 1) % n][0], ring[(i + 1) % n][1]
            if (y1 > lat) != (y2 > lat) and lon < (x2 - x1) * (lat - y1) / (y2 - y1 + 1e-15) + x1: inside = not inside
    return inside

def png_classes(b):
    a = np.array(Image.open(io.BytesIO(b)).convert("RGBA"))
    cl = np.full(a.shape[:2], -1, np.int8); rgb = a[..., :3].astype(int); op = a[..., 3] > 0
    for i, c in enumerate(COL):
        cl[op & (rgb[..., 0] == c[0]) & (rgb[..., 1] == c[1]) & (rgb[..., 2] == c[2])] = i
    odd = int((op & (cl < 0)).sum())
    if odd:   # 見知らぬ色（気象庁が色を変えたとき）: いちばん近い色にする
        d = np.stack([((rgb - np.array(c)) ** 2).sum(-1) for c in COL], -1); k = d.argmin(-1)
        m = op & (cl < 0); cl[m] = k[m]
    return cl, odd

def fetch_time(ts, mem="none"):
    """解析の時刻 ts（YYYYMMDDHHMM00）の一式: 1・3・24時間の 1km の階級と 5km の最大。読めなかったら None"""
    out = {"ts": ts, "odd": 0}
    for el in ("rasrf", "rasrf03h", "rasrf24h"):
        st, b = get(f"{D}{ts}/{mem}/{ts}/surf/{el}_point/data.geojson")
        if st != 200 or not b:
            if el == "rasrf": return None, f"{el}_point {st}"
            out[el] = None; out[el + "_m"] = None; continue   # 3・24時間がないときは、制約に使わないだけ
        try: fs = json.loads(b)["features"]
        except Exception as e: return None, f"{el}_point json {e}"
        cm = np.zeros((CJ1 - CJ0 + 1, CI1 - CI0 + 1), np.float32)
        for f in fs:
            m = re.search(r"_(\d+)_(\d+)$", f["properties"]["id"])
            i, j = int(m.group(1)), int(m.group(2))
            if CI0 <= i <= CI1 and CJ0 <= j <= CJ1: cm[j - CJ0, i - CI0] = float(f["properties"]["value"])
        grid = np.full((NY, NX), -1, np.int8); tiles = {}
        def one(xy):
            x, y = xy
            st, b = get(f"{D}{ts}/{mem}/{ts}/surf/{el}/10/{x}/{y}.png")
            return xy, st, b
        with ThreadPoolExecutor(6) as ex: res = list(ex.map(one, TL))
        for (xy, st, b) in res:
            if st == 200 and b[:4] == b"\x89PNG":
                cl, odd = png_classes(b); tiles[xy] = cl; out["odd"] += odd
            elif st == 404: tiles[xy] = None   # 雨のないタイル（気象庁は作らない）
            else:
                if el == "rasrf": return None, f"{el} tile {xy} {st}"
                tiles = None; break
        if tiles is None: out[el] = None; out[el + "_m"] = None; continue
        for k, xy in enumerate(TL):
            t = tiles[xy]
            if t is None: continue
            m = TIX == k
            grid[m] = t[PYS[m], PXS[m]]
        # つじつまの確かめ: 5km 四方に雨の値があるのに、その中の 1km 格子がどれも「雨なし」→ 画像がまだ出そろっていない（あとでもう一度）
        bad = 0
        for jj in range(cm.shape[0]):
            for ii in range(cm.shape[1]):
                if cm[jj, ii] <= 0: continue
                q0, p0 = (CJ0 + jj) * 6 - Q0, (CI0 + ii) * 5 - P0
                if q0 < 0 or p0 < 0 or q0 + 6 > NY or p0 + 5 > NX: continue
                if (grid[q0:q0 + 6, p0:p0 + 5] < 0).all(): bad += 1
        if bad:
            if el == "rasrf": return None, f"{el} tiles not ready ({bad} cells)"
            out[el] = None; out[el + "_m"] = None; continue
        out[el] = grid; out[el + "_m"] = cm
    # 解析のない所（レーダーの止まった所など）: 福井の陸にかかっていたら、この時間は「欠け」とする
    st, b = get(f"{D}{ts}/{mem}/{ts}/surf/rasrf_nd/data.geojson")
    out["nd"] = False
    if st == 200 and b:
        try:
            for f in json.loads(b).get("features", []):
                g = f.get("geometry") or {}
                polys = [g["coordinates"]] if g.get("type") == "Polygon" else g.get("coordinates", []) if g.get("type") == "MultiPolygon" else []
                for poly in polys:
                    if any(_in_rings(lat, lon, poly) for lat, lon in ND_PTS): out["nd"] = True
        except Exception as e:
            log("  nd parse", ts, e)
    return out, ""

def save_state(d, rec):
    kw = {"c1": rec["rasrf"], "m1": rec["rasrf_m"], "nd": np.array(bool(rec.get("nd")))}
    if rec.get("rasrf03h") is not None: kw.update(c3=rec["rasrf03h"], m3=rec["rasrf03h_m"])
    if rec.get("rasrf24h") is not None: kw.update(c24=rec["rasrf24h"], m24=rec["rasrf24h_m"])
    np.savez_compressed(os.path.join(d, rec["ts"] + ".npz"), **kw)
def load_state(d, ts):
    f = os.path.join(d, ts + ".npz")
    if not os.path.exists(f): return None
    try: z = np.load(f)
    except Exception: return None
    return {k: z[k] for k in z.files}

def mesh_max(cm): return cm[CELLJ, CELLI]

def estimate(recs):
    """recs: 古い→新しい順の1時間ごとの記録（None は欠け）。返り値: 1時間ごとの 1km 格子の雨量（mm）の配列 [H, NY, NX] と、欠けの印"""
    H = len(recs); X = np.zeros((H, NY, NX), np.float32); LO = np.zeros_like(X); HI = np.zeros_like(X)
    N0 = np.zeros((H, NY, NX), np.float32); miss = np.zeros(H, bool)
    for h, r in enumerate(recs):
        if r is None: miss[h] = True; continue
        if "nd" in r and bool(r["nd"]): miss[h] = True   # 解析のない時間（値は使うが、欠けとして数える）
        c = r["c1"]; M = mesh_max(r["m1"])
        lo = np.zeros((NY, NX), np.float32); hi = np.zeros((NY, NX), np.float32); x = np.zeros((NY, NX), np.float32)
        for k, (a, b) in enumerate(B1):
            m = c == k
            if not m.any(): continue
            if k == 0: x[m] = C0; lo[m] = C0; hi[m] = C0; continue
            bb = np.maximum(M[m], a) if b is None else np.full(m.sum(), b, np.float32)
            top = np.minimum(bb, np.where(M[m] > 0, M[m], bb))
            lo[m] = np.minimum(a, top); hi[m] = top; x[m] = np.minimum(a + FRAC * (bb - a), top)
        N0[h] = c == 0
        adj = c >= 1
        for _ in range(SMOOTH):
            pad = np.pad(x, 1, mode="edge"); s = sum(pad[1 + dy:1 + dy + NY, 1 + dx:1 + dx + NX] for dy in (-1, 0, 1) for dx in (-1, 0, 1)) / 9
            x = np.where(adj, np.clip(s, np.minimum(lo, x), hi), x)
        X[h] = x; LO[h] = np.where(adj, np.minimum(lo, x), x); HI[h] = np.where(adj, hi, x)
    # 3時間・24時間の制約（気象庁の積算＝1mm以上の時間の値の和＋0.4×1mm未満の時間の数）
    cons = []
    for w, key, BB in ((3, "3", B3), (24, "24", B24)):
        for h in range(w - 1, H):
            r = recs[h]
            if r is None or miss[h - w + 1:h + 1].any() or ("c" + key) not in r: continue
            c = r["c" + key]; M = mesh_max(r["m" + key])
            lo = np.zeros((NY, NX), np.float32); hi = np.full((NY, NX), 1e9, np.float32)
            for k, (a, b) in enumerate(BB):
                m = c == k; lo[m] = a; hi[m] = 1e9 if b is None else b
            hi[c == -1] = 1.0
            hi = np.minimum(hi, np.where(M > 0, M + 0.05, hi))
            n0 = N0[h - w + 1:h + 1].sum(0); x0 = (X[h - w + 1:h + 1] * N0[h - w + 1:h + 1]).sum(0)
            cons.append((h - w + 1, h, lo - 0.4 * n0 + x0, hi - 0.4 * n0 + x0))
    for _ in range(6):
        for (i0, i1, lo, hi) in cons:
            seg = X[i0:i1 + 1]; s = seg.sum(0)
            up = s < lo - 1e-4
            if up.any():
                room = (HI[i0:i1 + 1] - seg).sum(0); f = np.where(up & (room > 0), np.minimum(1, (lo - s) / np.maximum(room, 1e-9)), 0)
                seg += (HI[i0:i1 + 1] - seg) * f
            dn = s > hi + 1e-4
            if dn.any():
                room = (seg - LO[i0:i1 + 1]).sum(0); f = np.where(dn & (room > 0), np.minimum(1, (s - hi) / np.maximum(room, 1e-9)), 0)
                seg -= (seg - LO[i0:i1 + 1]) * f
            X[i0:i1 + 1] = seg
    return X, miss

def b64u16(a):
    a = np.clip(np.round(a * 10), 0, 65534).astype("<u2")
    return base64.b64encode(a.tobytes()).decode()

def mesh_of(lat, lon): return int((lon - 118) / 0.0125), int((48 - lat) * 120)

def check_amedas(T, X, miss, out_dir):
    """アメダス（福井とまわり）と答え合わせ。12・24・48時間"""
    st, b = get(A + "const/amedastable.json")
    if st != 200: return None
    table = json.loads(b)
    def deg(v): return v[0] + v[1] / 60
    sts = {}
    for k, v in table.items():
        if str(v.get("elems", "0"))[1:2] != "1": continue
        lat, lon = deg(v["lat"]), deg(v["lon"]); p, q = mesh_of(lat, lon)
        if P0 <= p < P1 and Q0 <= q < Q1: sts[k] = (v["kjName"], q - Q0, p - P0, lat, lon)
    tj = T + timedelta(hours=9)
    need = {}
    for code in sts:
        for back in (0, 3, 6, 9, 24):
            t = tj - timedelta(hours=back); bl = t.strftime("%Y%m%d_") + "%02d" % (t.hour // 3 * 3)
            need[(code, bl)] = None
    def one(k):
        code, bl = k; s, bb = get(f"{A}data/point/{code}/{bl}.json")
        return k, (json.loads(bb) if s == 200 else None)
    with ThreadPoolExecutor(8) as ex:
        for k, d in ex.map(one, list(need)): need[k] = d
    def rec(code, back, key):
        t = tj - timedelta(hours=back); bl = t.strftime("%Y%m%d_") + "%02d" % (t.hour // 3 * 3)
        d = need.get((code, bl)); r = d and d.get(t.strftime("%Y%m%d%H%M00"))
        v = r and r.get(key)
        return v[0] if v and v[0] is not None and v[1] == 0 else None
    rows = []
    for code, (nm, qi, pi, lat, lon) in sts.items():
        g12 = [rec(code, b, "precipitation3h") for b in (0, 3, 6, 9)]
        g = {12: None if None in g12 else sum(g12), 24: rec(code, 0, "precipitation24h")}
        g48 = [rec(code, 0, "precipitation24h"), rec(code, 24, "precipitation24h")]
        g[48] = None if None in g48 else sum(g48)
        e = {W: (None if miss[-W:].any() else float(X[-W:, qi, pi].sum())) for W in WINS}
        rows.append({"code": code, "name": nm, "lat": round(lat, 4), "lon": round(lon, 4), **{f"g{W}": g[W] for W in WINS}, **{f"e{W}": None if e[W] is None else round(e[W], 1) for W in WINS}})
    summ = {}
    for W in WINS:
        pr = [(r[f"e{W}"], r[f"g{W}"]) for r in rows if r[f"e{W}"] is not None and r[f"g{W}"] is not None and (r[f"g{W}"] >= 1 or r[f"e{W}"] >= 1)]
        if not pr: summ[str(W)] = {"n": 0}; continue
        e = np.array([p[0] for p in pr]); g = np.array([p[1] for p in pr]); ae = np.abs(e - g)
        summ[str(W)] = {"n": len(pr), "bias": round(float((e - g).mean()), 2), "mae": round(float(ae.mean()), 2), "ratio": round(float(e.sum() / max(g.sum(), 1e-9)), 3),
                        "med_rel": round(float(np.median(ae / np.maximum(g, 1))), 3), "ok20": round(float(((ae <= 2) | (ae / np.maximum(g, 1) <= 0.2)).mean()), 3)}
    return {"t": T.strftime("%Y-%m-%dT%H:%M:%SZ"), "summary": summ, "rows": rows}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--budget", type=float, default=900); ap.add_argument("--max-new", type=int, default=60)
    ap.add_argument("--at", default="")   # 試験用: この時刻で作る（YYYYMMDDHHMM00）
    a = ap.parse_args()
    t0 = time.time(); os.makedirs(a.state, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    st, b = get(D + "targetTimes.json")
    if st != 200: log("targetTimes", st); sys.exit(1)
    R = json.loads(b)
    an = sorted(e["validtime"] for e in R if e["basetime"] == e["validtime"] and (e.get("member") or "none") == "none" and "rasrf24h" in (e.get("elements") or []))
    if not an: log("no analysis"); sys.exit(1)
    T = datetime.strptime(a.at or an[-1], "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)
    series = [(T - timedelta(hours=k)).strftime("%Y%m%d%H%M00") for k in range(KEEP_H - 1, -1, -1)]   # 古い→新しい
    have = set(f[:-4] for f in os.listdir(a.state) if f.endswith(".npz"))
    todo = [ts for ts in reversed(series) if ts not in have][:a.max_new]   # 新しい時刻から
    log("latest", an[-1], "use", T.strftime("%Y%m%d%H%M"), "need", len([s for s in series if s not in have]), "fetch now", len(todo))
    failed = []
    ff = os.path.join(a.state, "fail.json")
    try: fails = json.load(open(ff))
    except Exception: fails = {}
    newest = series[-1]
    for ts in todo:
        if time.time() - t0 > a.budget: log("budget over"); break
        if fails.get(ts, 0) >= 3 and ts != newest: continue   # 3回読めなかった時間は、欠けのまま（毎回は読みに行かない）
        rec, why = fetch_time(ts)
        if rec is None:
            failed.append((ts, why)); log("  fail", ts, why)
            if "not ready" not in why: fails[ts] = fails.get(ts, 0) + 1
            continue
        save_state(a.state, rec); fails.pop(ts, None)
        if rec["odd"]: log("  odd colors", ts, rec["odd"])
        if rec.get("nd"): log("  no-data area over Fukui", ts)
    lim0 = (T - timedelta(hours=KEEP_H + 2)).strftime("%Y%m%d%H%M00")
    json.dump({k: v for k, v in fails.items() if k >= lim0}, open(ff, "w"))
    # 古い記録を消す（ほかの時刻の並び（:00 と :30）も、それぞれ KEEP_H 時間）
    lim = (T - timedelta(hours=KEEP_H + 2)).strftime("%Y%m%d%H%M00")
    for f in os.listdir(a.state):
        if f.endswith(".npz") and f[:-4] < lim: os.remove(os.path.join(a.state, f))
    recs = [load_state(a.state, ts) for ts in series]
    X, miss = estimate(recs)
    out = {"v": 1, "t": T.strftime("%Y-%m-%dT%H:%M:%SZ"), "made": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "src": "気象庁ホームページ（解析雨量）を加工して作成", "grid": {"p0": P0, "q0": Q0, "nx": NX, "ny": NY, "lon0": 118, "lat0": 48, "dlon": 0.0125, "dlat": 1 / 120},
           "w": list(WINS), "unit": 0.1, "a": {}, "miss": {}}
    for W in WINS:
        out["a"][str(W)] = b64u16(X[-W:].sum(0)); out["miss"][str(W)] = int(miss[-W:].sum())
    # 答え合わせ（3時間ごとの正時）
    chk = None
    if T.minute == 0 and T.hour % 3 == 0:
        try: chk = check_amedas(T, X, miss, a.out)
        except Exception as e: log("check failed", e)
        if chk:
            out["chk"] = chk["summary"]
            with open(os.path.join(a.out, "check.jsonl"), "a") as f: f.write(json.dumps(chk, ensure_ascii=False) + "\n")
            lines = open(os.path.join(a.out, "check.jsonl")).read().splitlines()[-3000:]
            open(os.path.join(a.out, "check.jsonl"), "w").write("\n".join(lines) + "\n")
            log("check", json.dumps(chk["summary"], ensure_ascii=False))
    json.dump(out, open(os.path.join(a.out, "acc.json"), "w"), ensure_ascii=False, separators=(",", ":"))
    open(os.path.join(a.out, "t.txt"), "w").write(out["t"] + "\n")
    log("done", out["t"], {W: out["miss"][str(W)] for W in WINS}, "max mm", {W: round(float(X[-W:].sum(0).max()), 1) for W in WINS}, "failed", len(failed), round(time.time() - t0), "s")

if __name__ == "__main__":
    main()
