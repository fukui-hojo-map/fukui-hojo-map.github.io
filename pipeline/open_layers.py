"""Earth Engine を使わない地図の層: 土地被覆（ESA WorldCover 2021）と標高（Copernicus DEM 30m）。AWS の公開データ（鍵なし）。

Earth Engine の `ESA/WorldCover/v200`・`COPERNICUS/DEM/GLO30` と同じデータを、同じ格子（pixels.grid_for の EPSG:3857）に写す。
  WorldCover: 最も近い画素（Earth Engine の computePixels と同じ）。3度四方のタイル（1/12000 度 ≒ 10m）
  DEM       : 両隣の画素からの直線補間（Earth Engine の resample("bilinear") と同じ考え）。1度四方のタイル（1秒 ≒ 30m）
wild.py（獣害の手がかり: 森＝10）・near.py（住宅地の近さ: 市街地＝50）・wet.py（地形の予備の標高）が使う口:
  OpenWild().forest(g)、OpenNear().built(g)、OpenDEM().dem(g)   （GEEWild・GEENear・GEEWet.dem と同じ形）
"""
import math
import numpy as np

WC = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{t}_Map.tif"
DEM = "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{t}_DEM/Copernicus_DSM_COG_10_{t}_DEM.tif"
R = 6378137.0


def lonlat_of(g):
    """格子 g（EPSG:3857）の画素の中心の経度・緯度（H×W）"""
    xs = g["x0"] + (np.arange(g["w"]) + 0.5) * g["px"]; ys = g["y1"] - (np.arange(g["h"]) + 0.5) * g["px"]
    lon = np.degrees(xs / R); lat = np.degrees(2 * np.arctan(np.exp(ys / R)) - math.pi / 2)
    return np.meshgrid(lon, lat)


def _read(url, c0, r0, c1, r1):
    """COG の窓（列 c0〜c1、行 r0〜r1、端を含まない）。ファイルの外は 0"""
    import rasterio
    from rasterio.windows import Window
    from s2open import _gdal_env
    for k in range(5):
        try:
            with _gdal_env(), rasterio.open(url + (f"?retry={k}" if k else "")) as s:
                out = np.zeros((r1 - r0, c1 - c0), s.dtypes[0])
                ca, ra, cb, rb = max(c0, 0), max(r0, 0), min(c1, s.width), min(r1, s.height)
                if cb > ca and rb > ra:
                    out[ra - r0:rb - r0, ca - c0:cb - c0] = s.read(1, window=Window(ca, ra, cb - ca, rb - ra))
                return out
        except Exception:
            if k == 4:
                raise
            import time
            time.sleep(3 * 2 ** k)


def worldcover(g):
    """格子 g の WorldCover の区分（10=樹林, 40=耕地, 50=市街地, …。0=データなし）。最も近い画素"""
    lon, lat = lonlat_of(g); out = np.zeros(lon.shape, np.uint8); n = 12000
    for la in range(int(math.floor(lat.min() / 3)) * 3, int(math.floor(lat.max() / 3)) * 3 + 1, 3):
        for lo in range(int(math.floor(lon.min() / 3)) * 3, int(math.floor(lon.max() / 3)) * 3 + 1, 3):
            m = (lat >= la) & (lat < la + 3) & (lon >= lo) & (lon < lo + 3)
            if not m.any():
                continue
            c = np.floor((lon[m] - lo) * n).astype(np.int64); r = np.floor((la + 3 - lat[m]) * n).astype(np.int64)
            t = f"{'N' if la >= 0 else 'S'}{abs(la):02d}{'E' if lo >= 0 else 'W'}{abs(lo):03d}"
            a = _read(WC.format(t=t), c.min(), r.min(), c.max() + 1, r.max() + 1)
            out[m] = a[r - r.min(), c - c.min()]
    return out


def dem(g):
    """格子 g の標高（m, float32。海・データなしは NaN）。直線補間（タイルの境目は、そのタイルの端の値で補う）"""
    from scipy import ndimage
    lon, lat = lonlat_of(g); out = np.full(lon.shape, np.nan, np.float32); n = 3600
    for la in range(int(math.floor(lat.min())), int(math.floor(lat.max())) + 1):
        for lo in range(int(math.floor(lon.min())), int(math.floor(lon.max())) + 1):
            m = (lat >= la) & (lat < la + 1) & (lon >= lo) & (lon < lo + 1)
            if not m.any():
                continue
            fc = (lon[m] - lo) * n - 0.5; fr = (la + 1 - lat[m]) * n - 0.5           # 画素の中心からの位置
            c0, r0 = int(math.floor(fc.min())), int(math.floor(fr.min())); c1, r1 = int(math.floor(fc.max())) + 2, int(math.floor(fr.max())) + 2
            t = f"{'N' if la >= 0 else 'S'}{abs(la):02d}_00_{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00"
            try:
                a = _read(DEM.format(t=t), max(c0, 0), max(r0, 0), min(c1, n), min(r1, n)).astype(np.float32)
            except Exception:
                continue                                                   # 海だけのタイルはない
            a = np.pad(a, ((max(0, -r0), max(0, r1 - n)), (max(0, -c0), max(0, c1 - n))), mode="edge")
            out[m] = ndimage.map_coordinates(a, [fr - r0, fc - c0], order=1, mode="nearest")
    return out


class OpenWild:
    """wild.GEEWild と同じ口: 森（WorldCover 10）"""
    def forest(self, g):
        return (worldcover(g) == 10).astype(np.uint8)


class OpenNear:
    """near.GEENear と同じ口: 市街地・建物（WorldCover 50）"""
    def built(self, g):
        return (worldcover(g) == 50).astype(np.uint8)


class OpenDEM:
    """wet.GEEWet.dem と同じ口"""
    def dem(self, g):
        return dem(g)
