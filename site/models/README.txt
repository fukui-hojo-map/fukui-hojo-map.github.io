AR画面の「モノを配置する」で置ける 3D モデル（開発者用）

■ モノを増やすとき
  1. この models/ フォルダに、.glb ファイルを入れる（できれば 10 MB 以下・三角形 10万個以下。大きいものは軽くする）
  2. site/index.html の AR_OBJECTS に、1行足す:
       {id:'名前(英数字)', name:'ボタンに出す名前', url:'models/ファイル名.glb', fit:['z',実物の長さ(m)], credit:'作者の表示（必要なら）'}
     fit は「モデルのどの軸（x・y・z）の長さを、実物の何mにするか」。置いたとき、モデルの前（−Z 側）が、置いた人のほうを向く（モデルの前が −Z でないときは、yaw0 で回す）
  3. ライセンスを確認する（CC BY は作者名の表示が必要・NC は商用利用ができない）。表示は credit に書くと、置いているあいだ画面に出る

■ 入っているもの
  lamborghini.glb
    "2021 Lamborghini Countach LPI 800-4" by Lexyc16
    https://sketchfab.com/3d-models/2021-lamborghini-countach-lpi-800-4-d76b94884432422b966d1a7f8815afb5
    CC BY-NC 4.0（作者の表示が必要・商用利用は不可）。ライセンスの原文: lamborghini.LICENSE.txt
    加工: 三角形を約26万個から約10万個に減らし、窓の透過を半透明に、テクスチャを 1024px の WebP にした
