function p_myfunc_drawFigure(mode, timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers)
     if mode == "trajectory"
        p_myfunc_drawTrajectoryQuiver(timeseries, U_sym)
     elseif mode == "transition_vec"
        p_myfunc_drawTransitionVecQuiver(timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval)
     elseif mode == "basins"
        p_myfunc_drawBasins(timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers)
     end
end


function p_myfunc_drawTrajectoryQuiver(timeseries, U_sym)

    figure;
    syms x1 x2
    x = linspace(-0.8,0.8);
    y = linspace(-0.7,0.7);
    [x,y] = meshgrid(x,y);
    U = matlabFunction(U_sym, 'Vars', {x1, x2});
    z = U(x,y);
    contour(x, y, z, 50);
    title('Surface Plot');
    % Show the colorbar for reference
    colorbar;
    hold on;

    x1_ = timeseries(:,1);
    x2_ = timeseries(:,2);
    
    % plot time series data
    plot(x1_, x2_, '-o', 'DisplayName', '時系列データ');
    
    % plot quiver
    for t = 1:length(x1_)-1
        quiver(x1_(t), x2_(t), x1_(t+1)-x1_(t), x2_(t+1)-x2_(t), 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 0.5);
    end
    
    % plot settings
    xlabel('x1');
    ylabel('x2');
    xlim([-0.8 0.8])
    ylim([-0.7 0.7])
    grid on;
    hold off;

end

function p_myfunc_drawTransitionVecQuiver(timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval)

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

    syms x1 x2

    figure;
    x = linspace(-0.8,0.8);
    y = linspace(-0.7,0.7);
    [x,y] = meshgrid(x,y);
    U = matlabFunction(U_sym, 'Vars', {x1, x2});
    z = U(x,y);
    contour(x, y, z, 50);
    title('Surface Plot');
    % Show the colorbar for reference
    colorbar;
    hold on;

    %[x_center, y_center] = meshgrid(1:num_x1_cells, 1:num_x2_cells); % セルの中心点のx座標とy座標を生成
    %x_center = x_center(:); % ベクトルに変換
    %y_center = y_center(:); % ベクトルに変換
    %x_center = x1_grid(x_center) + epsilon/2; % セルの中心点のx座標を計算
    %y_center = x2_grid(y_center) + epsilon/2; 
    
    % 計算した平均を元に、それぞれのセルごとに矢印を書くという操作で可視化する
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            if numel(average_vecs{i,j}) > 0
                %quiver(x_center((i-1)*num_x2_cells + j), y_center((i-1)*num_x2_cells + j), average_vecs{i,j}(1)/3, average_vecs{i,j}(2)/3, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 3)
                quiver(average_vec_start_points{i, j}(1), average_vec_start_points{i, j}(2), average_vecs{i,j}(1)/3, average_vecs{i,j}(2)/3, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 3)
                hold on;
            end
        end
    end
    xlabel('x1');
    ylabel('x2');
    xlim([-0.8 0.8])
    ylim([-0.7 0.7])
    for i = 1:num_x1_cells
        for j = 1:num_x2_cells
            rectangle('Position', [x1_grid(i), x2_grid(j), gridded_interval, gridded_interval], 'EdgeColor', 'k'); % グリッドをプロット
        end
    end
    grid on;
end

function p_myfunc_drawBasins(timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers)
    [m, n] = size(basins);
    unique_basins = unique(basins);
    num_basins = length(unique_basins);
    
    p_myfunc_drawTransitionVecQuiver(timeseries, average_vecs, average_vec_start_points, U_sym, gridded_interval)
    hold on;
    
    % 各 basin ごとにプロットとテキストの設定
    for i = 1:num_basins
        basin_id = unique_basins(i);
        if basin_id == 0
            continue; % basin_idが0の場合はスキップ
        end
        basin_indices = find(basins == basin_id);
        for j = 1:length(basin_indices)
            [row, col] = ind2sub([m, n], basin_indices(j));
            center = cell_centers{row, col};
            
            % プロット
            plot(center(1), center(2), 'o', 'MarkerSize', 0.000001, 'MarkerEdgeColor', 'k');
            
            % テキストの追加
            text(center(1), center(2), num2str(basin_id), 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
        end
    end
    
    hold off;
    xlabel('X');
    ylabel('Y');
    title('Basins Classification');
    colorbar; % カラーバーの追加
  end