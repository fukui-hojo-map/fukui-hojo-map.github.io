#!/bin/bash
# ある置き場（ディレクトリ）の中身を、同じ名前のブランチに置く。変わったときだけ、履歴をためずに1つのコミットで置きかえる
#   bash pipeline/branch_push.sh <ディレクトリ> <ブランチ> <コミットの文>
# 環境変数: GITHUB_TOKEN, GITHUB_REPOSITORY（RADAR_REMOTE_URL は試験用に置き先を変える）。置き場からファイルが減るときは置かない（取り出しに失敗している）
set -e
dir="$1"; branch="$2"; msg="$3"
cd "$dir"
[ -d .git ] || git init -q
git config user.name github-actions
git config user.email github-actions@users.noreply.github.com
url="${RADAR_REMOTE_URL:-https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git}"
git remote set-url origin "$url" 2>/dev/null || git remote add origin "$url"
git config --unset-all http.https://github.com/.extraheader || true
git add -A
if git rev-parse -q --verify HEAD >/dev/null; then
  if git diff --cached --quiet; then echo "$branch: 変わりなし"; exit 0; fi
  if [ -n "$(git diff --cached --diff-filter=D --name-only)" ]; then echo "$branch: ファイルが減っているので置きません"; exit 1; fi
  git checkout -q --orphan "${branch}-next"      # 履歴をためない（毎回1つのコミットにする）
  git add -A
fi
git commit -qm "$msg"
git push -q -f origin "HEAD:${branch}"
echo "$branch: 置きました"
