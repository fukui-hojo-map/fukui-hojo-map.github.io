#!/bin/bash
# レーダーのデータを radar-data ブランチに置く（pipeline/radar_store.py の置き場 radar/ を、変わったときだけ1つのコミットで置きかえる）
# weekly-ndvi から呼ぶ。失敗したら終了コードが 0 以外（呼び出し側は、置けなかったときはレーダーを公開サイトに残す）
set -e
python pipeline/radar_store.py save --data site/data --store radar
bash pipeline/branch_push.sh radar radar-data "radar data"
