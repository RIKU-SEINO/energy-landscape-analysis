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

  % energy_values を新しい定義に基づいて更新する
  for i = 1:num_x1_cells
    for j = 1:num_x2_cells
      if ~isempty(cell_vecs{i, j})
        % 指定された式に基づいて新しいエネルギー値を計算する
        exp_term = 0;
        for k = 1:num_x1_cells
          for l = 1:num_x2_cells
            if ~isempty(cell_vecs{k, l})
              exp_term = exp_term + exp(-2 * energy_values{k, l} / sigma^2);
            end
          end
        end
        energy_values{i, j} = - (sigma^2) / 2 * log(probability_values{i, j}) + log(exp_term);
      end
    end
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
  colorbar;

  U = 200*(0.2*X.^4 + 0.4*Y.^4 - 0.1*X.^2 - 0.1*Y.^2);

  hold on; % 現在のプロットに追加
  surf(X, Y, U, 'FaceAlpha', 0.5); % 半透明で U をプロット
  xlabel('x1');
  ylabel('x2');
  zlabel('U');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  title('出現確率の分布と U の3次元プロット');
  colorbar;
  hold off;