function [probability_values, energy_values] = p_myfunc_computeEnergyFromProbabilityDist(cell_vecs, sigma)
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
        energy_values{i, j} = 0;
      end
    end
  end

  % 3次元バー グラフを作成してプロットする
  figure;
  dataArray = cellfun(@double, probability_values); 
  x = linspace(-0.8, 0.8, size(dataArray, 2));
  y = linspace(-0.7, 0.7, size(dataArray, 1));
  [X, Y] = meshgrid(x,y);
  surf(X, Y, dataArray);
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('P(x1,x2)');
  title('出現確率の分布');

  figure;
  dataArray = cellfun(@double, energy_values); 
  x = linspace(-0.8, 0.8, size(dataArray, 2));
  y = linspace(-0.7, 0.7, size(dataArray, 1));
  [X, Y] = meshgrid(x,y);
  surf(X, Y, dataArray);
  xlabel('x1');
  ylabel('x2');
  xlim([-0.8 0.8])
  ylim([-0.7 0.7])
  zlabel('E(x1,x2)');
  title('エネルギーの分布');


  