function result = p_myfunc_cosSim(average_vecs, average_vec_start_points)

    [num_x1_cells, num_x2_cells] = size(average_vecs);
    cos_sims = cell(num_x1_cells, num_x2_cells);
    cnt = 0;
    cos_sim_sum = 0;
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            if numel(average_vecs{i,j}) > 0
                average_vec = average_vecs{i,j};
                average_vec_start_point = average_vec_start_points{i,j};
                gradient_at_average_vec_start_point = F(0, average_vec_start_point);%Fはtに依存しないように設定しているので、t=0としている、別になんでもいいけど
                cos_sim = dot(average_vec, gradient_at_average_vec_start_point) / (norm(average_vec) * norm(gradient_at_average_vec_start_point));
                cos_sims{i,j} = cos_sim;
                cos_sim_sum = cos_sim_sum + cos_sim;
                cnt = cnt + 1;
            end
        end
    end

    result = cos_sim_sum / cnt;
