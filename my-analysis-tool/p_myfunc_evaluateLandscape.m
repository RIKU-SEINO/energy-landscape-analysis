function energy_values_real = p_myfunc_evaluateLandscape(energy_values, energy_values_real, x1_grid, x2_grid)

  % セル配列を作成し、実際のエネルギー地形と構築したエネルギー地形の差異をセルごとに格納
  energy_values_diff = cell(size(energy_values));
  for i = 1:size(energy_values, 1)
      for j = 1:size(energy_values, 2)
          energy_values_diff{i, j} = energy_values_real{i, j} - energy_values{i, j};
      end
  end

  % 3次元バー グラフを作成してプロットする
  figure;
  dataArray = cellfun(@double, energy_values_diff); 
  [X, Y] = meshgrid(x1_grid, x2_grid);
  surf(X, Y, dataArray.');
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('P(x1,x2)');
  title('エネルギー差の分布');

  % カラーバーを追加
  colorbar;
