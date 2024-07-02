function [result1, result2] = p_myfunc_cosSim(average_vecs, average_vec_start_points, cell_vec_start_points, U_sym)
    % U_sym をシンボリック変数として定義
    syms x1 x2
    
    % U_sym を x1 と x2 で微分
    grad_U = gradient(U_sym, [x1, x2]);
    
    % grad_U を数値関数として定義
    grad_U_func = matlabFunction(grad_U, 'Vars', {x1, x2});

    % average_vecs のサイズを取得
    [num_x1_cells, num_x2_cells] = size(average_vecs);
    
    % コサイン類似度のセル配列を初期化
    cos_sims1 = cell(num_x1_cells, num_x2_cells);
    cos_sims2 = cell(num_x1_cells, num_x2_cells);
    
    % 計算用の変数を初期化
    cnt = 0;
    cos_sim1_sum = 0;
    cos_sim2_sum = 0;
    
    % average_vecs のループ
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            % average_vecs{i,j} が空でない場合
            if numel(average_vecs{i,j}) > 0
                average_vec = average_vecs{i,j};
                average_vec_start_point = average_vec_start_points{i,j};
                
                % average_vec_start_point での勾配を計算
                gradient_at_average_vec_start_point = -grad_U_func(average_vec_start_point(1), average_vec_start_point(2));

                % vec_start_pointsそれぞれにおける勾配を計算し、それを平均とる
                vec_start_points = cell_vec_start_points{i, j};
                for ii = 1:size(vec_start_points, 2)
                    vec_start_point = vec_start_points(:,ii);
                    average_gradient_at_each_vec_start_point = -grad_U_func(vec_start_point(1), vec_start_point(2)) / size(vec_start_points, 2);
                end
                
                % コサイン類似度を計算
                norm1 = norm(average_vec);
                norm2 = norm(gradient_at_average_vec_start_point);
                norm3 = norm(average_gradient_at_each_vec_start_point);
                if norm1 == 0 || norm2 == 0 || norm3 == 0
                    continue;
                end
                cos_sim1 = dot(average_vec, gradient_at_average_vec_start_point) / (norm1 * norm2);
                cos_sim2 = dot(average_vec, average_gradient_at_each_vec_start_point) / (norm1 * norm3);
                % 結果をセル配列に格納
                cos_sims1{i,j} = cos_sim1;
                cos_sims2{i,j} = cos_sim2;
                cos_sim1_sum = cos_sim1_sum + cos_sim1;
                cos_sim2_sum = cos_sim2_sum + cos_sim2;
                cnt = cnt + 1;
            end
        end
    end

    % 平均コサイン類似度を計算
    result1 = cos_sim1_sum / cnt;
    result2 = cos_sim2_sum / cnt;
