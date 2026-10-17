"""レーダー（Sentinel-1）のデータを、公開サイト（GitHub Pages）の外に置く。

公開サイトは約1GB の上限に近い。そのうちレーダーが約230MB（区画ごとの系列 ndvi.json の rp・rdates・rmask が約150MB、
過去の年の radar.json が約85MB）で、画面は読んでいない（rmask が古い形式のときだけ使う作りで、いまのデータでは使われない）。
そこで、レーダーだけ別のブランチ（radar-data）に置き、公開サイトには出さない。作物の見分け（crop.py）などの計算は、
これまでどおり公開前のデータ（site/data）に全部そろった形で行う。

  python pipeline/radar_store.py restore --data site/data --store radar
      別ブランチ（を取り出した radar/）の分を、site/data の ndvi.json などに戻す。更新の最初に1回
  python pipeline/radar_store.py save    --data site/data --store radar
      site/data のレーダー分を radar/ に書き出す。変わったファイルの数を出力する（0 なら別ブランチは書きかえない）
  python pipeline/radar_store.py strip   --data site_pub/data
      公開用の複製（cp -al site site_pub で作る。ハードリンク）から、レーダー分を取り除く。元の site/data は変えない

置き場の形（radar/）:
  cells/<id>/ndvi_r.json       {"v":2, "rmask", "rdates", "rp"}  ndvi.json から切り出した分
  cells/<id>/ndvi_next_r.json  同じ（方式が変わって計算し直し中の ndvi_next.json の分）
  cells/<id>/radar.json        過去の年（ndvi.json より前）のレーダー。そのまま
最初の1回は、公開中の ndvi.json にレーダーが入ったままなので、そのまま save すれば置き場ができる。
restore は、公開側にレーダーが残っているときはそれを優先する（置き場より新しい可能性がある）。
"""
import os, json, argparse
from common import read_json

KEYS = ("rmask", "rdates", "rp")
MAIN = ("ndvi.json", "ndvi_next.json")        # レーダーを含みうる系列ファイル
WHOLE = ("radar.json",)                       # そのまま置くファイル


def r_name(fname):
    return fname[:-5] + "_r.json"


def cell_dirs(root):
    base = os.path.join(root, "cells")
    if not os.path.isdir(base):
        return []
    return sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))


def dumps(obj):
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def write_text(path, text):
    """書く。ハードリンクで他と同じ実体になっていても、元を壊さないよう、いったん消してから書く"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.lexists(path):
        os.remove(path)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def has_radar(nd):
    return bool(nd) and any(k in nd for k in KEYS)


def piece(nd):
    """ndvi.json のうちレーダー分（置き場に置く形）"""
    out = {"v": 2}
    for k in KEYS:
        if k in nd:
            out[k] = nd[k]
    return out


def restore(data, store):
    """置き場 → site/data。公開側にレーダーが残っていれば、そちらを使う"""
    n = 0
    for cid in cell_dirs(store):
        sdir, ddir = os.path.join(store, "cells", cid), os.path.join(data, "cells", cid)
        if not os.path.isdir(ddir):
            continue
        for fname in MAIN:
            r = read_json(os.path.join(sdir, r_name(fname)))
            nd = read_json(os.path.join(ddir, fname))
            if not r or not nd or has_radar(nd):
                continue
            nd.update({k: r[k] for k in KEYS if k in r})
            write_text(os.path.join(ddir, fname), dumps(nd)); n += 1
        for fname in WHOLE:
            p = os.path.join(sdir, fname)
            if os.path.exists(p) and not os.path.exists(os.path.join(ddir, fname)):
                with open(p, encoding="utf-8") as f:
                    write_text(os.path.join(ddir, fname), f.read()); n += 1
    return n


def save(data, store):
    """site/data → 置き場。内容が変わったファイルだけ書く。書いた数を返す"""
    n = 0
    for cid in cell_dirs(data):
        ddir, sdir = os.path.join(data, "cells", cid), os.path.join(store, "cells", cid)
        for fname in MAIN:
            nd = read_json(os.path.join(ddir, fname))
            if not has_radar(nd):
                continue
            text = dumps(piece(nd)); p = os.path.join(sdir, r_name(fname))
            if not os.path.exists(p) or open(p, encoding="utf-8").read() != text:
                write_text(p, text); n += 1
        for fname in WHOLE:
            p = os.path.join(ddir, fname)
            if not os.path.exists(p):
                continue
            with open(p, encoding="utf-8") as f:
                text = f.read()
            q = os.path.join(sdir, fname)
            if not os.path.exists(q) or open(q, encoding="utf-8").read() != text:
                write_text(q, text); n += 1
    return n


def strip(data):
    """公開用の複製から、レーダー分を取り除く。取り除いたファイル数を返す"""
    n = 0
    for cid in cell_dirs(data):
        ddir = os.path.join(data, "cells", cid)
        for fname in MAIN:
            p = os.path.join(ddir, fname)
            nd = read_json(p)
            if not has_radar(nd):
                continue
            for k in KEYS:
                nd.pop(k, None)
            write_text(p, dumps(nd)); n += 1
        for fname in WHOLE:
            p = os.path.join(ddir, fname)
            if os.path.exists(p):
                os.remove(p); n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="レーダーのデータを公開サイトの外に置く")
    ap.add_argument("cmd", choices=["restore", "save", "strip"])
    ap.add_argument("--data", required=True, help="データの置き場（site/data）")
    ap.add_argument("--store", help="レーダー用の置き場（別ブランチを取り出した所）")
    a = ap.parse_args()
    if a.cmd in ("restore", "save") and not a.store:
        ap.error("--store が要ります")
    n = {"restore": lambda: restore(a.data, a.store), "save": lambda: save(a.data, a.store), "strip": lambda: strip(a.data)}[a.cmd]()
    print(f"{a.cmd}: {n} ファイル")


if __name__ == "__main__":
    main()
