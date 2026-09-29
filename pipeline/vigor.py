"""生育傾向マップの田ごとの値（例年の育ち）を作る: cells/<id>/vig.json と vig_stats.json（Earth Engine は使わない）。

福井は早生（ハナエチゼン）から晩生、直播まで作期がばらばらで、同じ日の NDVI をくらべると「早く植えた・早く育つ品種」の
田ほど高く出る（6月ごろの差の大半は植えた時期の差だった）。そこで田ごと・稲の年ごとに
  頭打ち vmax = 6/14〜9/19 の晴れた日（区画の内側の8割以上）の NDVI の最大（その年いちばん茂った時。出穂のころ）
  立ち上がり t50 = 頭打ちより前で NDVI が 0.5 を上に越えた日（植えた時期の目安）
  刈り取り th = 頭打ちから10日以上あとで 0.45 未満に落ちた最初の日と、その前の観測の中間（熟期の目安。見えなければなし）
を出し、近く（3×3 セル ≒ 7×5.5 km）の、t50 と th が同じ6日の枠に入る田（＝同じころに植えて同じころに刈った田）の
平均とくらべる（その田を除く。10枚に足りなければ枠を ±1 に広げ、さらに 5×5 セルに広げる）。
小さい田は縁（畦・道）が画素に混ざって低く出るので、区画の内側の面積ごとの差（大きい田＝0）を先に差し引く（年ごとにデータから求める）。
田の値 = 年ごとの差の平均 × k·σb²/(k·σb²+σw²)（k = 年数。1〜2年だけの田が極端に出ないよう、年ごとのぶれの分だけ 0 に寄せる）。
  2022〜2026年の県全体で: 年どうしの一致 0.58（田の値の約半分は毎年くり返す差）、熟期（刈り取り日）との相関 0.08、
  植えた時期との相関 −0.12（今までの 6月ごろ・8月ごろ の値は −0.50）、100m 以内の田どうしの相関 0.49。

  python pipeline/vigor.py [--data DIR] [--out DIR] [--prev DIR] [--today YYYY-MM-DD] [--workers 4]

--prev DIR: 前回の出力。ndvi.json（いつもの系列は24か月で古い日を落とす）や hist.json から消えた年は、前回の vig.json に
            残した頭打ち・立ち上がり・刈り取りを使う（くらべる相手は毎回計算し直す）。
今年は 9/20 から入れる（晩生もいちばん茂る時期を過ぎる。まだ刈り取りが見えない田は、見えない田どうしでくらべる。刈り取りが見えるたびに値が少し動く）。
cells/<id>/vig.json: {"v":1, "years":[年...], "p":{pid:[田の値×1000|null, [年, 差×1000|null, 頭打ち×1000, くらべた田の数, t50, th|0], ...]}}
vig_stats.json: {"v":1, "years", "sb", "sw", "q":[県全体の田の値×1000 の 0〜100% 点], "n", "size":{年:[内側の面積の区切りごとの差×1000]}}
"""
import argparse, datetime, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), 'site', 'data')
JST = datetime.timezone(datetime.timedelta(hours=9))
VERSION = 1
CLEAR_MIN = 80                 # 区画の内側の晴れた割合（%）
PEAK_WIN = (165, 262)          # 頭打ちをさがす日（通日。6/14〜9/19）
HARV_MAX, HARV_AFTER = 0.45, 10
UP = 0.5
BIN = 6                        # t50・th の枠（日）
NEED = 10                      # くらべる田の数の下限
LEVELS = ((1, 0, 0), (1, 1, 1), (2, 1, 1))   # (セルの半径, t50 の枠の広げ, th の枠の広げ)
SIZE_EDGES = [0, 100, 200, 300, 400, 500, 650, 800, 1000, 1250, 1500, 2000, 2500, 3000, 4000, 6000, 1e12]
SIZE_BASE = 3000               # この面積（㎡、区画の内側）以上の田を 0 にそろえる
MIN_YEARS = 2


def doy(s):
    d = datetime.date(int(s[:4]), int(s[5:7]), int(s[8:10]))
    return d.year, (d - datetime.date(d.year, 1, 1)).days + 1


def ring_area(r, lat0):
    kx = 111320 * math.cos(math.radians(lat0)); ky = 110540; s = 0.0
    for i in range(len(r) - 1):
        s += (r[i][0] * kx) * (r[i + 1][1] * ky) - (r[i + 1][0] * kx) * (r[i][1] * ky)
    return abs(s) / 2


def geom_area(gm):
    polys = gm['coordinates'] if gm['type'] == 'MultiPolygon' else [gm['coordinates']]
    lat0 = sum(p[1] for p in polys[0][0]) / len(polys[0][0])
    return sum(ring_area(pl[0], lat0) - sum(ring_area(h, lat0) for h in pl[1:]) for pl in polys)


def read(path):
    try:
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None


def season_feats(t, v):
    """晴れた日の系列（通日・NDVI、日付順）→ (頭打ち, t50, th) か None"""
    m = (t >= PEAK_WIN[0]) & (t <= PEAK_WIN[1])
    if m.sum() < 2:
        return None
    i = np.flatnonzero(m)[np.argmax(v[m])]; tmax, vmax = t[i], v[i]
    th = None
    after = np.flatnonzero((t > tmax + HARV_AFTER) & (v < HARV_MAX))
    if len(after):
        j = after[0]; th = (t[j] + t[j - 1]) / 2 if j > 0 else float(t[j])
    k = t < (th if th is not None else 1e9); tt, vv = t[k], v[k]
    b = np.flatnonzero((tt < tmax) & (vv < UP))
    if not len(b) or b[-1] + 1 >= len(tt):
        return None
    j = b[-1]
    t50 = tt[j] + (UP - vv[j]) / (vv[j + 1] - vv[j]) * (tt[j + 1] - tt[j]) if vv[j + 1] > vv[j] else tt[j + 1]
    return float(vmax), float(t50), th


def load_cell(args):
    """1セル: 稲の年の田ごとの (年, 頭打ち, t50, th) と、区画の内側の面積"""
    cid, data, prev, years = args
    cd = os.path.join(data, 'cells', cid)
    cs = read(os.path.join(cd, 'crop_s.json'))
    if not cs:
        return cid, {}, {}, 'crop_s.json なし'
    rice = {}
    for pid, a in cs['p'].items():
        ys = {y for y, q in zip(cs['years'], a) if q and q[0] == '稲' and y in years}
        if ys:
            rice[pid] = ys
    if not rice:
        return cid, {}, {}, ''
    ser = {}
    for fn in ('hist.json', 'post.json', 'ndvi.json'):              # post.json: 2022〜2024年の 10/1〜11/30（10月刈りの田の刈り取り）
        h = read(os.path.join(cd, fn))
        if not h or not h.get('p'):
            continue
        dts = [doy(s) for s in h['dates']]
        for pid, a in h['p'].items():
            if pid not in rice:
                continue
            for k in range(0, len(a), 3):
                y, d = dts[a[k]]
                if y in rice[pid] and a[k + 2] >= CLEAR_MIN and 121 <= d <= 304:
                    ser.setdefault((pid, y), {})[d] = a[k + 1] / 1000     # 同じ日が2つのファイルにあれば後（ndvi.json）を使う
    feats = {}
    for (pid, y), dv in ser.items():
        t = np.array(sorted(dv), float); v = np.array([dv[d] for d in sorted(dv)], float)
        f = season_feats(t, v)
        if f:
            feats[(pid, y)] = f
    pv = read(os.path.join(prev, 'cells', cid, 'vig.json')) if prev else None
    if pv and pv.get('v') == VERSION:                          # 元の系列から消えた年は前回の値を使う
        for pid, a in pv['p'].items():
            for q in a[1:]:
                y = q[0]
                if pid in rice and y in rice[pid] and (pid, y) not in feats and (pid, y) not in ser:
                    feats[(pid, y)] = (q[2] / 1000, float(q[4]), float(q[5]) if q[5] else None)
    area = {}
    inner = read(os.path.join(cd, 'inner.geojson'))
    for f in (inner or {}).get('features', []):
        pid = f['properties'].get('pid')
        if pid in rice and f.get('geometry'):
            area[pid] = geom_area(f['geometry'])
    return cid, feats, area, ''


def _nsum(A, R, dt, dh):
    P = np.pad(A, ((R, R), (R, R), (dt, dt), (dh, dh))); out = np.zeros_like(A); nx, ny, nt, nh = A.shape
    for a in range(2 * R + 1):
        for b in range(2 * R + 1):
            for c in range(2 * dt + 1):
                for d in range(2 * dh + 1):
                    out += P[a:a + nx, b:b + ny, c:c + nt, d:d + nh]
    return out


def peer_resid(val, yr, ix, iy, t50, th):
    """近くの、t50・th が同じ枠の田の平均との差（その田を除く。値は年ごとに上下1%で切ってから平均）→ (差, くらべた田の数)"""
    res = np.full(len(val), np.nan); cnt = np.zeros(len(val), int)
    ix = ix - ix.min(); iy = iy - iy.min()
    bt = np.floor(t50 / BIN).astype(int); bt -= bt.min()
    has = ~np.isnan(th); bh = np.zeros(len(val), int)
    bh[has] = np.floor(th[has] / BIN).astype(int); bh[has] -= bh[has].min() if has.any() else 0
    nh = bh.max() + 1; bh[~has] = nh + 2                        # 刈り取りが見えない田は、見えた田と混ぜない
    shp = (ix.max() + 1, iy.max() + 1, bt.max() + 1, nh + 3)
    for y in np.unique(yr):
        k = yr == y
        lo, hi = np.percentile(val[k], [1, 99]); vv = np.clip(val[k], lo, hi)
        S = np.zeros(shp); N = np.zeros(shp); key = (ix[k], iy[k], bt[k], bh[k])
        np.add.at(S, key, vv); np.add.at(N, key, 1)
        r = np.full(k.sum(), np.nan); c = np.zeros(k.sum(), int)
        for R, dt, dh in LEVELS:
            todo = np.isnan(r)
            if not todo.any():
                break
            Sn = _nsum(S, R, dt, dh)[key] - vv; Nn = _nsum(N, R, dt, dh)[key] - 1
            ok = todo & (Nn >= NEED)
            r[ok] = val[k][ok] - Sn[ok] / Nn[ok]; c[ok] = Nn[ok]
        res[k] = r; cnt[k] = c
    return res, cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=DATA); ap.add_argument('--out', default=DATA); ap.add_argument('--prev', default='')
    ap.add_argument('--today', default=''); ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--from-year', type=int, default=2022)
    a = ap.parse_args()
    t0 = time.time()
    today = datetime.date.fromisoformat(a.today) if a.today else datetime.datetime.now(JST).date()
    last = today.year if today >= datetime.date(today.year, 9, 20) else today.year - 1
    years = list(range(a.from_year, last + 1))
    idx = read(os.path.join(a.data, 'index.json'))
    cells = [c['id'] for c in idx['cells']]
    rows = []; area = {}; errs = []
    with ProcessPoolExecutor(a.workers) as ex:
        for cid, feats, ar, err in ex.map(load_cell, [(c, a.data, a.prev, set(years)) for c in cells], chunksize=4):
            if err:
                errs.append(f'{cid}: {err}')
            for (pid, y), (vmax, t50, th) in feats.items():
                rows.append((cid, pid, y, vmax, t50, np.nan if th is None else th))
            area.update(ar)
    if not rows:
        print('稲の年のデータがありません'); return 1
    cid_ = np.array([r[0] for r in rows]); pid_ = np.array([r[1] for r in rows]); yr = np.array([r[2] for r in rows])
    vmax = np.array([r[3] for r in rows]); t50 = np.array([r[4] for r in rows]); th = np.array([r[5] for r in rows])
    ix = np.array([int(c.split('_')[0]) for c in cid_]); iy = np.array([int(c.split('_')[1]) for c in cid_])
    ai = np.array([area.get(p, np.nan) for p in pid_])
    # 1回目: 面積を考えずにくらべる → 内側の面積ごとの差の中央値（年ごと。大きい田＝0）→ 差し引いて、くらべ直す
    r0, _ = peer_resid(vmax, yr, ix, iy, t50, th)
    sb_ = np.digitize(np.nan_to_num(ai), SIZE_EDGES) - 1
    adj = np.zeros(len(rows)); size = {}
    for y in np.unique(yr):
        k = yr == y; med = np.full(len(SIZE_EDGES) - 1, np.nan)
        for b in range(len(SIZE_EDGES) - 1):
            q = k & (sb_ == b) & ~np.isnan(r0)
            if q.sum() >= 50:
                med[b] = np.median(r0[q])
        base = np.nanmean(med[[b for b in range(len(med)) if SIZE_EDGES[b] >= SIZE_BASE]])
        med = np.where(np.isnan(med), 0, med - (0 if np.isnan(base) else base))
        adj[k] = med[sb_[k]]; size[int(y)] = [int(round(x * 1000)) for x in med]
    res, cnt = peer_resid(vmax - adj, yr, ix, iy, t50, th)
    # 年ごとのぶれ（σw）と田どうしの差（σb）→ 縮め
    ok = ~np.isnan(res)
    byp = {}
    for p, r in zip(pid_[ok], res[ok]):
        byp.setdefault(p, []).append(r)
    multi = [np.array(v) for v in byp.values() if len(v) >= 2]
    sw2 = sum(((v - v.mean()) ** 2).sum() for v in multi) / max(1, sum(len(v) - 1 for v in multi))
    fm = np.array([v.mean() for v in multi]); fn = np.array([len(v) for v in multi])
    sb2 = max(1e-5, fm.var() - np.mean(sw2 / fn))
    score = {p: float(np.mean(v)) * len(v) * sb2 / (len(v) * sb2 + sw2) for p, v in byp.items() if len(v) >= MIN_YEARS}
    # 書き出し
    per = {}
    for i in range(len(rows)):
        per.setdefault(cid_[i], {}).setdefault(pid_[i], []).append(
            [int(yr[i]), None if np.isnan(res[i]) else int(round(res[i] * 1000)), int(round(vmax[i] * 1000)), int(cnt[i]),
             int(round(t50[i])), 0 if np.isnan(th[i]) else int(round(th[i]))])
    for c in cells:
        p = per.get(c, {})
        out = {'v': VERSION, 'years': years,
               'p': {pid: [None if pid not in score else int(round(score[pid] * 1000))] + sorted(v) for pid, v in p.items()}}
        d = os.path.join(a.out, 'cells', c); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'vig.json'), 'w', encoding='utf-8') as f:
            json.dump(out, f, ensure_ascii=False, separators=(',', ':'))
    sc = np.array(list(score.values()))
    stats = {'v': VERSION, 'updated': datetime.datetime.now(JST).isoformat(timespec='seconds'), 'years': years,
             'sb': round(math.sqrt(sb2), 4), 'sw': round(math.sqrt(sw2), 4), 'n': len(sc),
             'q': [int(round(x * 1000)) for x in np.percentile(sc, np.arange(101))] if len(sc) else [],
             'size_edges': SIZE_EDGES[:-1], 'size': size}
    with open(os.path.join(a.out, 'vig_stats.json'), 'w', encoding='utf-8') as f:
        json.dump(stats, f, ensure_ascii=False, separators=(',', ':'))
    print(f'田・年 {len(rows)}（くらべられた {int(ok.sum())}）、田の値 {len(sc)} 枚、年 {years}、σb {stats["sb"]} σw {stats["sw"]}、'
          f'{time.time() - t0:.0f} 秒' + (f'、読めなかったセル {len(errs)}' if errs else ''))
    for e in errs[:10]:
        print('  ', e)
    return 2 if len(errs) > 0.05 * len(cells) else 0


if __name__ == '__main__':
    sys.exit(main())
