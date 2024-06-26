function [average_vecs, average_vec_start_points, variance_vecs, average_vec_lengths, counts] = p_myfunc_statsForEachCell(cell_vecs, cell_vec_start_points)

    [num_x1_cells, num_x2_cells] = size(cell_vecs);
    average_vecs = cell(num_x1_cells, num_x2_cells); % 平均化された変化ベクトルを格納するセル配列を初期化
    average_vec_lengths = cell(num_x1_cells, num_x2_cells); %平均化された変化ベクトルのノルムを格納するセル配列を初期化
    variance_vecs = cell(num_x1_cells, num_x2_cells); % 変化ベクトルの分散を格納するセル配列を初期化
    average_vec_start_points = cell(num_x1_cells, num_x2_cells); % 変化ベクトルの始点を平均化して得られた座標を格納するセル配列を初期化
    counts = cell(num_x1_cells, num_x2_cells); % 変化ベクトルの数を格納するセル配列を初期化
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            % セル内の変化ベクトルを取得
            cell_vec = cell_vecs{i, j};
            cell_vec_start_point = cell_vec_start_points{i, j};
            % セル内に変化ベクトルが存在する場合
            if ~isempty(cell_vec)
                % 変化ベクトルの始点の平均を計算
                average_vec_start_points{i, j} = mean(cell_vec_start_point, 2);
                % 平均化された変化ベクトルを計算
                average_vecs{i, j} = mean(cell_vec,2);
                % 分散を計算
                variance_vecs{i, j} = cov(cell_vec');
                % 平均化された変化ベクトルのノルムを計算
                average_vec_lengths{i, j} = norm(average_vec_lengths{i, j}, 2);
                % カウントを計算
                counts{i, j} = length(cell_vec);
            end
        end
    end