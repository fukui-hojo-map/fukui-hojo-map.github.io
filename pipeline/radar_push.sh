#!/bin/bash
# レーダーのデータを radar-data ブランチに置く（pipeline/radar_store.py の置き場 radar/ を、変わったときだけ1つのコミットで置きかえる）
# weekly-ndvi から呼ぶ。環境変数: GITHUB_TOKEN, GITHUB_REPOSITORY。失敗したら終了コードが 0 以外（呼び出し側は、置けなかったときはレーダーを公開サイトに残す）
set -e
python pipeline/radar_store.py save --data site/data --store radar
cd radar
[ -d .git ] || git init -q
git config user.name github-actions
git config user.email github-actions@users.noreply.github.com
url="${RADAR_REMOTE_URL:-https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git}"   # RADAR_REMOTE_URL は試験用
git remote set-url origin "$url" 2>/dev/null || git remote add origin "$url"
git config --unset-all http.https://github.com/.extraheader || true
git add -A
if git rev-parse -q --verify HEAD >/dev/null; then
  if git diff --cached --quiet; then echo "radar-data: 変わりなし"; exit 0; fi
  # 置き場からファイルが消えることはない（save は足すだけ）。消えるなら取り出しに失敗している。置きかえない
  if [ -n "$(git diff --cached --diff-filter=D --name-only)" ]; then echo "radar-data: ファイルが減っているので置きません"; exit 1; fi
  git checkout -q --orphan radar-next      # 履歴をためない（毎回1つのコミットにする）
  git add -A
fi
git commit -qm "radar data"
git push -q -f origin HEAD:radar-data
echo "radar-data: 置きました"
