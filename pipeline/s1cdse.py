"""Earth Engine を使わないレーダー（Sentinel-1 の VH）。Copernicus Data Space Ecosystem（CDSE）の Sentinel Hub から取る。

Earth Engine の COPERNICUS/S1_GRD（update_ndvi.GEEBackend の src="s1"）と同じ考えにそろえる:
  - IW・VV+VH（DV）・高解像度（HIGH）。σ0（楕円体で正規化。backCoeff SIGMA0_ELLIPSOID）、地形の補正（orthorectify、Copernicus DEM 30m）
  - 同じ日（UTC）の画像は1枚に重ねる（新しい方が上）
  - 区画の値は、内側の画素（10m）の VH（線形）の中央値。数は画素の数。値は小数7桁（update_ndvi が dB にする）
  - 生育傾向の R（trend.py）は、期間（5/1〜6/10）の画素ごとの VH（dB）の中央値 → 区画の中央値
格子は s2open と同じ（UTM、10m、9.6km のブロック）。区画の重みも s2open の Parcels・weights をそのまま使う。

鍵: 環境変数 CDSE_CLIENT_ID・CDSE_CLIENT_SECRET（Sentinel Hub のダッシュボードで作る OAuth クライアント。GitHub のシークレットに置く。チャット・コード・ログに出さない）
使う量: Sentinel Hub の「処理の単位（PU）」と「回数」に、無料の月ごとの上限がある。1ブロック・1日を1回で取り（1バンド・16bit）、結果は置き場に残して、同じものは2回取らない。
  取った回数と、Sentinel Hub が返す使った PU（ヘッダ）を stats_n に数える（試しの実行で、月の使う量を見積もる）

  python pipeline/s1cdse.py dates --start 2025-05-01 --end 2025-05-15        （鍵があるとき: 撮った日の一覧だけ。PU は使わない）
"""
import os, io, json, time, threading, urllib.request, urllib.error, urllib.parse, datetime as dt
import numpy as np

TOKEN_URL = os.environ.get("CDSE_TOKEN_URL") or "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
SH_URL = os.environ.get("CDSE_SH_URL") or "https://sh.dataspace.copernicus.eu"
COLLECTION = "sentinel-1-grd"
VERSION = "c1"            # 取り方の版。変えると置き場の結果を取り直す
DB0, DBK = 50.0, 1000.0   # VH の dB を 16bit で受け取る: 値 = round((dB + 50) × 1000)、0 = データなし（0.001 dB きざみ）
EVALSCRIPT = """//VERSION=3
function setup(){ return {input:[{bands:["VH","dataMask"]}], output:{bands:1, sampleType:"UINT16"}}; }
function evaluatePixel(s){ if(!s.dataMask||!(s.VH>0)) return [0]; const db=10*Math.log(s.VH)/Math.LN10; return [Math.max(1,Math.min(65535,Math.round((db+%s)*%s)))]; }
""" % (DB0, DBK)


class CDSE:
    """CDSE の Sentinel Hub（Catalog・Process）を呼ぶ。トークンは期限の1分前に取り直す。429・5xx は待ってやり直す"""
    def __init__(self, client_id=None, client_secret=None, opener=None):
        self.cid = client_id or os.environ.get("CDSE_CLIENT_ID") or ""
        self.sec = client_secret or os.environ.get("CDSE_CLIENT_SECRET") or ""
        self.open = opener or urllib.request.urlopen
        self._tok = None; self._exp = 0.0; self._lock = threading.Lock()
        self.stats_n = {"process": 0, "catalog": 0, "pu": 0.0, "retry": 0}

    def ok(self):
        return bool(self.cid and self.sec)

    def token(self):
        with self._lock:
            if self._tok and time.time() < self._exp - 60:
                return self._tok
            if not self.ok():
                raise RuntimeError("CDSE の鍵（CDSE_CLIENT_ID・CDSE_CLIENT_SECRET）がありません")
            body = urllib.parse.urlencode({"grant_type": "client_credentials", "client_id": self.cid, "client_secret": self.sec}).encode()
            req = urllib.request.Request(TOKEN_URL, data=body, headers={"content-type": "application/x-www-form-urlencoded"})
            with self.open(req, timeout=60) as r:
                j = json.loads(r.read())
            self._tok = j["access_token"]; self._exp = time.time() + float(j.get("expires_in") or 600)
            return self._tok

    def _post(self, url, obj, accept="application/json", tries=6):
        for k in range(tries):
            req = urllib.request.Request(url, data=json.dumps(obj).encode(), method="POST",
                                         headers={"authorization": "Bearer " + self.token(), "content-type": "application/json", "accept": accept})
            try:
                with self.open(req, timeout=180) as r:
                    pu = r.headers.get("x-processingunits-spent") if hasattr(r, "headers") and r.headers else None
                    if pu:
                        try:
                            self.stats_n["pu"] += float(pu)
                        except ValueError:
                            pass
                    return r.read()
            except urllib.error.HTTPError as e:
                code = e.code
                if code == 401 and k == 0:                      # トークンが切れていた
                    with self._lock:
                        self._tok = None
                    continue
                if k == tries - 1 or (code < 500 and code != 429):
                    msg = e.read()[:300].decode("utf-8", "replace") if hasattr(e, "read") else ""
                    raise RuntimeError(f"CDSE {code}: {msg}")
                ra = e.headers.get("retry-after") if e.headers else None
                self.stats_n["retry"] += 1
                time.sleep(min(60.0, float(ra)) if ra and ra.replace(".", "", 1).isdigit() else 2 * 2 ** k)
            except (urllib.error.URLError, TimeoutError) as e:
                if k == tries - 1:
                    raise
                self.stats_n["retry"] += 1; time.sleep(2 * 2 ** k)

    def search(self, bbox, start, end):
        """撮った画像（IW・VH あり）の一覧 [{id, t（ミリ秒）, date（UTC）, geom}]。bbox は経緯度、end は含まない"""
        out = []; body = {"bbox": list(bbox), "datetime": f"{start}T00:00:00Z/{end}T00:00:00Z", "collections": [COLLECTION], "limit": 100}
        for _ in range(100):
            j = json.loads(self._post(SH_URL + "/api/v1/catalog/1.0.0/search", body, accept="application/geo+json")); self.stats_n["catalog"] += 1   # application/json だと 406（2026-11-05 の試し）
            for f in j.get("features", []):
                p = f.get("properties", {})
                mode = p.get("sar:instrument_mode"); pol = p.get("sar:polarizations") or []
                if (mode and mode != "IW") or (pol and "VH" not in pol):
                    continue
                t = dt.datetime.fromisoformat(p["datetime"].replace("Z", "+00:00"))
                out.append({"id": f["id"], "t": int(t.timestamp() * 1000), "date": t.astimezone(dt.timezone.utc).strftime("%Y-%m-%d"), "geom": f.get("geometry")})
            nxt = (j.get("context") or {}).get("next")
            if nxt is None:
                break
            body = dict(body, next=nxt)
        return sorted({x["id"]: x for x in out}.values(), key=lambda x: x["t"])

    def block(self, date, x0, y1, size_m, px, epsg):
        """その日（UTC）の、UTM の四角（左の x0・上の y1・一辺 size_m）の VH（dB, float32。データなしは NaN）。px×px 画素"""
        nd = (dt.date.fromisoformat(date) + dt.timedelta(days=1)).isoformat()
        req = {"input": {"bounds": {"bbox": [x0, y1 - size_m, x0 + size_m, y1], "properties": {"crs": f"http://www.opengis.net/def/crs/EPSG/0/{epsg}"}},
                         "data": [{"type": COLLECTION,
                                   "dataFilter": {"timeRange": {"from": f"{date}T00:00:00Z", "to": f"{nd}T00:00:00Z"}, "acquisitionMode": "IW",
                                                  "polarization": "DV", "resolution": "HIGH", "mosaickingOrder": "mostRecent"},
                                   "processing": {"backCoeff": "SIGMA0_ELLIPSOID", "orthorectify": True, "demInstance": "COPERNICUS_30"}}]},
               "output": {"width": px, "height": px, "responses": [{"identifier": "default", "format": {"type": "image/tiff"}}]},
               "evalscript": EVALSCRIPT}
        raw = self._post(SH_URL + "/api/v1/process", req, accept="image/tiff"); self.stats_n["process"] += 1
        return decode_tiff(raw)


def decode_tiff(raw):
    """Process API の GeoTIFF（UINT16、0 = なし）→ dB（float32、NaN = なし）"""
    from rasterio.io import MemoryFile
    with MemoryFile(raw) as mf, mf.open() as ds:
        a = ds.read(1)
    v = a.astype(np.float32) / DBK - DB0; v[a == 0] = np.nan
    return v


class OpenS1:
    """update_ndvi.GEEBackend（src="s1"）と同じ口の、レーダーの部分。s2open.OpenS2（区画・重み・格子）と組んで使う"""
    def __init__(self, eng, api=None, cache=None):
        self.eng = eng; self.api = api or CDSE()
        self.cache = cache or os.path.join(eng.cache, "s1cdse")
        self._items = {}; self._lock = threading.Lock(); self._klocks = {}
        self.stats_n = {"blocks": 0, "cache_hit": 0}

    def _scenes(self, bbox, start, end):
        """bbox にかかる画像。県全体（eng.bbox）で1回だけ調べ、セルごとには写っている範囲で選ぶ（セルの数だけ問い合わせない）"""
        area = getattr(self.eng, "bbox", None)
        key = (start, end) if area else (tuple(round(x, 4) for x in bbox), start, end)
        with self._lock:
            if key not in self._items:
                self._items[key] = self.api.search(area or bbox, start, end)
            items = self._items[key]
            if not area:
                return items
            from shapely.geometry import shape, box
            b = box(*bbox); out = []
            for it in items:
                if it.get("geom") and "_shape" not in it:
                    it["_shape"] = shape(it["geom"])
                if "_shape" not in it or it["_shape"].intersects(b):
                    out.append(it)
            return out

    def images(self, bbox, start, end):
        return [{"id": it["id"], "t": it["t"], "geom": it["geom"]} for it in self._scenes(bbox, start, end)]

    def list_dates(self, bbox, start, end):
        return sorted({it["date"] for it in self._scenes(bbox, start, end)})

    def block_db(self, date, bkey):
        """ブロックの芯（960×960、10m）の、その日の VH（dB、NaN = なし）。置き場にあればそれ"""
        from s2open import BLOCK_PX, RES
        path = os.path.join(self.cache, VERSION, date, f"{bkey[0]}_{bkey[1]}.npz")
        with self._lock:
            kl = self._klocks.setdefault((date, bkey), threading.Lock())
        with kl:
            if os.path.exists(path):
                self.stats_n["cache_hit"] += 1
                z = np.load(path); a = z["q"]
                v = a.astype(np.float32) / DBK - DB0; v[a == 0] = np.nan
                return v
            L = BLOCK_PX * RES
            v = self.api.block(date, bkey[0] * L, (bkey[1] + 1) * L, L, BLOCK_PX, self.eng.epsg)
            q = np.where(np.isnan(v), 0, np.clip(np.round((v + DB0) * DBK), 1, 65535)).astype(np.uint16)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = path + f".{os.getpid()}.{threading.get_ident()}.tmp.npz"
            np.savez_compressed(tmp, q=q); os.replace(tmp, path)
            self.stats_n["blocks"] += 1
            return v

    def _parcel_vals(self, gis, img_of):
        """区画（添字）ごとの [(値, 重み, マスの数)…] を、ブロックごとの画像（芯、NaN = なし）から集める → (G, V, W, C)。重み・数は Earth Engine と同じ集め方（s2open.EE_DEG）"""
        from s2open import WIN_PX, HALO_PX
        P = self.eng.parcels; want = set(int(x) for x in gis); G, V, W, C = [], [], [], []
        for bk in P.blocks_of(gis):
            img = img_of(bk)
            if img is None:
                continue
            g, pix, w, wc = self.eng.weights(bk)
            m = np.fromiter((x in want for x in g.tolist()), bool, len(g))
            if not m.any():
                continue
            r = pix[m] // WIN_PX - HALO_PX; c = pix[m] % WIN_PX - HALO_PX
            v = img[r, c]; ok = ~np.isnan(v)
            G.append(g[m][ok]); V.append(v[ok]); W.append(w[m][ok]); C.append(wc[m][ok])
        if not G:
            return np.zeros(0, np.int64), np.zeros(0), np.zeros(0), np.zeros(0)
        return np.concatenate(G), np.concatenate(V), np.concatenate(W), np.concatenate(C)

    def stats(self, bbox, inner_fc_geojson, dates):
        """{日: {pid: [VH（線形）の中央値（小数7桁）, 画素の数]}}（GEEBackend.stats(src="s1") と同じ形）"""
        from trend import wmedian
        P = self.eng.parcels
        pids = [f["properties"]["pid"] for f in inner_fc_geojson["features"] if f["properties"]["pid"] in P.index]
        gis = [P.index[p] for p in pids]; out = {}
        have = set(self.list_dates(bbox, min(dates), (dt.date.fromisoformat(max(dates)) + dt.timedelta(days=1)).isoformat())) if dates else set()
        for d in dates:
            out[d] = {}
            if not gis or d not in have:
                continue
            G, V, W, C = self._parcel_vals(gis, lambda bk: self.block_db(d, bk))
            if not len(G):
                continue
            lin = 10 ** (V / 10)                                   # 中央値は線形でも dB でも同じ順番なので、どちらでとってもよい（線形で返す）
            med = wmedian(G, lin, W); cnt = np.bincount(G, weights=C, minlength=int(G.max()) + 1)   # 数は使えたマスの数（Earth Engine の count）
            for p, gi in zip(pids, gis):
                if gi in med:
                    out[d][p] = [round(float(med[gi]), 7), int(round(cnt[gi]))]
        return out

    def period_median_db(self, bkey, dates):
        """期間の日の、画素ごとの VH（dB）の中央値（生育傾向の R 用）。日がなければ None"""
        import warnings
        st = [self.block_db(d, bkey) for d in dates]
        if not st:
            return None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmedian(np.stack(st), axis=0).astype(np.float32)


def main():
    import argparse, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["dates"])
    ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--bbox", default="135.4,35.3,136.9,36.5")
    a = ap.parse_args()
    api = CDSE(); items = api.search([float(x) for x in a.bbox.split(",")], a.start, a.end)
    for it in items:
        print(it["date"], it["id"])
    print(len(items), "枚", api.stats_n)


if __name__ == "__main__":
    main()
