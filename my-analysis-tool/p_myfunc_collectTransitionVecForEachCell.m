function [cell_vecs, cell_vec_start_points] = p_myfunc_collectTransitionVecForEachCell(timeseries, vecs, gridded_interval)

    x1 = timeseries(:,1);
    x2 = timeseries(:,2);

    x1_ss_min = min(x1);
    x2_ss_min = min(x2);
    x1_ss_max = max(x1);
    x2_ss_max = max(x2);

    x1_grid = x1_ss_min:gridded_interval:x1_ss_max;
    x2_grid = x2_ss_min:gridded_interval:x2_ss_max;
    num_x1_cells = length(x1_grid);
    num_x2_cells = length(x2_grid);
    
    % cellごとに、変化ベクトルの集合vecs と 変化ベクトルの始点座標の集合vec_start_point を格納
    cell_vecs = cell(num_x1_cells, num_x2_cells);%各セル内に存在する変化ベクトルをセルごとに全てを格納、そのための初期化
    cell_vec_start_points = cell(num_x1_cells, num_x2_cells);% 各セル内に存在する変化ベクトルの始点をセルごとに全てを格納、そのための初期化
    for i = 1:length(vecs)
        x1_now = x1(i);
        x2_now = x2(i);
        x1_index = find(x1_now >= x1_grid, 1, 'last');
        x2_index = find(x2_now >= x2_grid, 1, 'last');
        if isempty(cell_vecs{x1_index, x2_index})
            cell_vecs{x1_index, x2_index} = [vecs(:,i)];
            cell_vec_start_points{x1_index, x2_index} = [x1_now;x2_now];
        else
            cell_vecs{x1_index, x2_index} = [cell_vecs{x1_index, x2_index}, vecs(:,i)];
            cell_vec_start_points{x1_index, x2_index} = [cell_vec_start_points{x1_index, x2_index}, [x1_now;x2_now]];
        end
    end