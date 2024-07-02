function energy_values_real = p_myfunc_constructLandscape(U_sym, interval_x1, interval_x2, energy_values)
  % シンボリック変数の定義
  syms x1 x2

  % 数値計算用の関数に変換
  U = matlabFunction(U_sym, 'Vars', [x1, x2]);

  % x1とx2の範囲を設定
  x1_range = linspace(-0.8, 0.8, interval_x1);
  x2_range = linspace(-0.7, 0.7, interval_x2);

  % メッシュグリッドを作成
  [X1, X2] = meshgrid(x1_range, x2_range);

  % ポテンシャル関数の値を計算
  U_vals = U(X1, X2);

  % セル配列を作成し、エネルギー値を格納
  energy_values_real = cell(size(U_vals));
  for i = 1:size(U_vals, 1)
      for j = 1:size(U_vals, 2)
          energy_values_real{i, j} = U_vals(i, j);
      end
  end

  % セル配列を作成し、実際のエネルギー地形と構築したエネルギー地形の差異をセルごとに格納
  energy_values_diff = cell(size(U_vals));
  for i = 1:size(U_vals, 1)
      for j = 1:size(U_vals, 2)
          energy_values_diff{i, j} = energy_values_real{i, j} - energy_values{i, j};
      end
  end


  % 3Dプロットを作成
  figure;
  surf(X1, X2, U_vals);

  % グラフのタイトルと軸ラベルを追加
  title('3D Plot of Potential Function');
  xlabel('x1');
  ylabel('x2');
  zlabel('U(x1, x2)');

  % 3次元バー グラフを作成してプロットする
  figure;
  dataArray = cellfun(@double, energy_values_diff); 
  x = linspace(-0.8, 0.8, size(dataArray, 2));
  y = linspace(-0.7, 0.7, size(dataArray, 1));
  [X, Y] = meshgrid(x,y);
  surf(X, Y, dataArray);
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('P(x1,x2)');
  title('エネルギー差の分布');

  % カラーバーを追加
  colorbar;
