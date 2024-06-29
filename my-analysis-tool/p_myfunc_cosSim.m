function result = p_myfunc_cosSim(average_vecs, average_vec_start_points, U_sym)
    % U_sym をシンボリック変数として定義
    syms x1 x2
    
    % U_sym を x1 と x2 で微分
    grad_U = gradient(U_sym, [x1, x2]);
    
    % grad_U を数値関数として定義
    grad_U_func = matlabFunction(grad_U, 'Vars', {x1, x2});

    % average_vecs のサイズを取得
    [num_x1_cells, num_x2_cells] = size(average_vecs);
    
    % コサイン類似度のセル配列を初期化
    cos_sims = cell(num_x1_cells, num_x2_cells);
    
    % 計算用の変数を初期化
    cnt = 0;
    cos_sim_sum = 0;
    
    % average_vecs のループ
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            % average_vecs{i,j} が空でない場合
            if numel(average_vecs{i,j}) > 0
                average_vec = average_vecs{i,j};
                average_vec_start_point = average_vec_start_points{i,j};
                
                % average_vec_start_point での勾配を計算
                gradient_at_average_vec_start_point = grad_U_func(average_vec_start_point(1), average_vec_start_point(2));
                
                % コサイン類似度を計算
                norm1 = norm(average_vec);
                norm2 = norm(gradient_at_average_vec_start_point);
                if norm1 == 0 || norm2 == 0
                    continue;
                end
                cos_sim = dot(average_vec, gradient_at_average_vec_start_point) / (norm1 * norm2);
                % 結果をセル配列に格納
                cos_sims{i,j} = cos_sim;
                cos_sim_sum = cos_sim_sum + cos_sim;
                cnt = cnt + 1;
            end
        end
    end

    % 平均コサイン類似度を計算
    result = cos_sim_sum / cnt;
