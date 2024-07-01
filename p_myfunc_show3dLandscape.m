function p_myfunc_show3dLandscape()
  % シンボリック変数の定義
  syms x1 x2

  % ポテンシャル関数の定義
  U_sym = 200 * (0.2 * x1^4 + 0.4 * x2^4 - 0.1 * x1^2 - 0.1 * x2^2);

  % 数値計算用の関数に変換
  U = matlabFunction(U_sym, 'Vars', [x1, x2]);

  % x1とx2の範囲を設定
  x1_range = linspace(-0.8, 0.8, 100);
  x2_range = linspace(-0.6, 0.6, 100);

  % メッシュグリッドを作成
  [X1, X2] = meshgrid(x1_range, x2_range);

  % ポテンシャル関数の値を計算
  U_vals = U(X1, X2);

  % 3Dプロットを作成
  figure;
  surf(X1, X2, U_vals);

  % グラフのタイトルと軸ラベルを追加
  title('3D Plot of Potential Function');
  xlabel('x1');
  ylabel('x2');
  zlabel('U(x1, x2)');

  % カラーバーを追加
  colorbar;
