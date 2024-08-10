function [probability_values, energy_values] = p_myfunc_computeEnergyFromProbabilityDist(cell_vecs, sigma, x1_grid, x2_grid)
  [num_x1_cells, num_x2_cells] = size(cell_vecs);
  probability_values = cell(num_x1_cells, num_x2_cells);
  energy_values = cell(num_x1_cells, num_x2_cells);

  total_num = 0;
  for i = 1:num_x1_cells
    for j = 1:num_x2_cells
      vecs = cell_vecs{i, j};
      if ~isempty(vecs)
        total_num = total_num + size(vecs, 2);
      end
    end
  end

  for i = 1:num_x1_cells
    for j = 1:num_x2_cells
      vecs = cell_vecs{i, j};
      if ~isempty(vecs)
        probability_values{i, j} = size(vecs, 2) / total_num;
        energy_values{i, j} = -(sigma.^2)*log(probability_values{i, j}) / 2;
      else
        probability_values{i, j} = 0;
        energy_values{i, j} = inf;
      end
    end
  end

  tol = 1e-6;  % 収束判定の閾値
  max_iter = 1000;  % 最大繰り返し回数
  iter = 0;
  delta = inf;  % 初期差分
  while delta > tol && iter < max_iter
    iter = iter + 1;
    old_energy_values = energy_values;  % エネルギー値を保存
    
    % エネルギーの更新
    C = log(sum(sum(exp(-2 * cellfun(@double, energy_values) / sigma^2))));
    for i = 1:num_x1_cells
      for j = 1:num_x2_cells
        if probability_values{i, j} > 0
          energy_values{i, j} = -((sigma.^2)/2) * (log(probability_values{i, j}) +  C);
        end
      end
    end
    
    % 収束判定
    delta = max(max(abs(cellfun(@double, energy_values) - cellfun(@double, old_energy_values))));
  end

  if iter >= max_iter
    warning('最大繰り返し回数に達しました。収束しなかった可能性があります。');
  end

  % 3次元バー グラフを作成してプロットする
  figure;
  dataArray = cellfun(@double, probability_values); 
  [X, Y] = meshgrid(x1_grid, x2_grid);
  surf(X, Y, dataArray.');
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('P(x1,x2)');
  title('出現確率の分布');

  figure;
  dataArray = cellfun(@double, energy_values); 
  [X, Y] = meshgrid(x1_grid, x2_grid);
  surf(X, Y, dataArray.');
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('E(x1,x2)');
  title('エネルギーの分布');