# ふくい圃場マップ

福井県内の圃場（筆ポリゴン）ごとの生育（NDVI）を、衛星（Sentinel-2）の画像から計算して地図で見せるウェブページです。

- 地図: https://fukui-hojo-map.github.io/
- つながるかの確認（行政端末など、通信先がしぼられた端末向け）: https://fukui-hojo-map.github.io/check.html

この置き場には、地図の画面（`site/`）と、公開データを作る自動処理の手順（`.github/workflows/`）だけを置いています。
**公開データを作る計算のプログラムは非公開**です（自動処理が、実行のときだけ非公開の置き場から取り出して使います）。

## 出典
- 衛星: Contains modified Copernicus Sentinel data。Sentinel-2 L2A・Sentinel-1 は Google Earth Engine（研究・非営利の登録）で計算（2026年11月25日から。指導者用のツールとして使うあいだ）。Landsat 8/9（NASA HLS）は 2026年10月までの分（Earth Engine で計算）
- 2026年11月14日〜25日は Earth Engine を使わずに計算した（Sentinel-2 L2A は AWS Open Data〔Element 84 の sentinel-cogs〕、Sentinel-1 は Copernicus Data Space Ecosystem）。その間の田の NDVI は、Earth Engine で計算し直したものに置きかえた
- 雲の判定: Cloud Score+（CC BY 4.0）。Earth Engine を使わない計算では、Cloud Score+ の判断を正解にして学習した判定器を使う。Pasquarella, V. J., Brown, C. F., Czerwinski, W., Rucklidge, W. J. (2023) "Comprehensive Quality Assessment of Optical Satellite Imagery Using Weakly Supervised Video Learning", CVPR Workshops 2023, doi:10.1109/CVPRW59228.2023.00206
- 区画: 農林水産省「筆ポリゴン」
- 地図: 国土地理院（地理院タイル〔航空写真・淡色地図〕・地名のベクトルタイル）。標高・高低差マップ: 国土地理院 標高タイル（基盤地図情報 数値標高モデル）を加工して作成
- 重ねる地図（洪水・土砂災害・津波の想定）: ハザードマップポータルサイト（国土交通省・国土地理院）
- 雨雲・雨量: 気象庁ホームページ（高解像度降水ナウキャスト・降水短時間予報・解析雨量）
- 土壌: 農研機構 日本土壌インベントリー（https://soil-inventory.rad.naro.go.jp/）CC BY 4.0
- 耕地・森: ESA WorldCover 2021（CC BY 4.0）。雪・雨: ERA5-Land（Copernicus Climate Change Service）

## 著作権
Copyright (c) 2026 山影祐也. All rights reserved.

この置き場にあるプログラム・画面・図・文章を、許可なく複製・改変・転載・再配布・転用することを禁じます（無断転用禁止）。くわしくは `LICENSE`。
地図に表示しているデータは、それぞれの提供元の利用条件に従います。
