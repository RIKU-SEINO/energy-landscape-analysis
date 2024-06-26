clear all
close all
%Create the following SDE composed of GBM and Mean Rverting drfit.
% dX(1) = 0.25*X(1)dt + 0.3*X(1)dW(1)
% dX(2) = 0.2*(0.1 -X(2))dt + 0 0.05*X(2)^0.5dW(2)

% define energy landscape function
figure(1)
x = linspace(-0.8,0.8);
y = linspace(-0.7,0.7);
[x,y] = meshgrid(x,y);
U = @(x,y)200*(0.2*x.^4+0.4*y.^4-0.1*x.^2-0.1*y.^2);
z = U(x,y);
contour(x, y, z, 50);
title('Surface Plot');
% Show the colorbar for reference
colorbar;



% U(x1,x2) = 200*(0.2*x1.^4 + 0.4*x2.^4 -0.1*x1.^2 -0.1*x2.^2);
% gradient of U = [(4*x1.^3)/5 - x1/5, (8*x2^3)/5 - x2/5]

% Define the system of equations
F = @(t,X) 200*[ - (4*X(1).^3)/5 + X(1)/5; -(8*X(2).^3)/5 + X(2)/5 ];
G = @(t,X) [ 0.8, 0; 0, 0.8];

% Define simulation parameters
startTime = 10;
%initialState = [1/2; -1/(2*sqrt(2));];  % Initial values for X(1) and X(2)
initialState = [0; 0];

% Create the sde object
obj = sde(F, G, 'StartTime', startTime, 'StartState', initialState);

% Simulate the system
nPeriods = 5000;
dt = 0.01;
stationaryStartIndex = 500;
[S,T] = simByEuler(obj, nPeriods, 'DeltaTime', dt);
x1 = S(:,1);
x2 = S(:,2);
x1_stationary = x1(stationaryStartIndex:end);
x2_stationary = x2(stationaryStartIndex:end);

hold on;

% plot time series data
plot(x1_stationary, x2_stationary, '-o', 'DisplayName', '時系列データ');

% plot quiver
for t = 1:length(x1_stationary)-1
    quiver(x1_stationary(t), x2_stationary(t), x1_stationary(t+1)-x1_stationary(t), x2_stationary(t+1)-x2_stationary(t), 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 0.5);
end

% plot settings
xlabel('x1');
ylabel('x2');
xlim([-0.8 0.8])
ylim([-0.7 0.7])
grid on;
hold off;

% output full .dat file
matrix = [x1; x2];
fid = fopen('fullTimeSeries.dat', 'w');
optionString = "%f\n";
for i = 1:length(x2)-1
    optionString = '%f ' + optionString;
end
fprintf(fid, optionString, matrix');
fclose(fid);

% output stationary .dat file
matrix = [x1_stationary; x2_stationary];
fid = fopen('stationaryTimeSeries.dat', 'w');
optionString = "%f\n";
for i = 1:length(x2_stationary)-1
    optionString = '%f ' + optionString;
end
fprintf(fid, optionString, matrix');
fclose(fid);

len = length(x1_stationary);
vecs = [];
for i = 1:len-1
    x_now = [x1_stationary(i);x2_stationary(i)];
    x_next = [x1_stationary(i+1);x2_stationary(i+1)];
    dx = x_next - x_now;
    vecs = [vecs, dx];
end

x1_ss_min = min(x1_stationary);
x2_ss_min = min(x2_stationary);
x1_ss_max = max(x1_stationary);
x2_ss_max = max(x2_stationary);
epsilon = 0.05;
x1_grid = x1_ss_min:epsilon:x1_ss_max;
x2_grid = x2_ss_min:epsilon:x2_ss_max;
num_x1_cells = length(x1_grid);
num_x2_cells = length(x2_grid);

% cellごとに、変化ベクトルの集合vecs と 変化ベクトルの始点座標の集合vec_start_point を格納
cell_vecs = cell(num_x1_cells, num_x2_cells);%各セル内に存在する変化ベクトルをセルごとに全てを格納、そのための初期化
cell_vec_start_points = cell(num_x1_cells, num_x2_cells);% 各セル内に存在する変化ベクトルの始点をセルごとに全てを格納、そのための初期化
for i = 1:length(vecs)
    x1_now = x1_stationary(i);
    x2_now = x2_stationary(i);
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

% ↑でセルごとに格納した情報を、セルごとに平均や分散を計算する
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

figure(2)
contour(x, y, z, 50);
title('Surface Plot');
% Show the colorbar for reference
colorbar;
hold on

%[x_center, y_center] = meshgrid(1:num_x1_cells, 1:num_x2_cells); % セルの中心点のx座標とy座標を生成
%x_center = x_center(:); % ベクトルに変換
%y_center = y_center(:); % ベクトルに変換
%x_center = x1_grid(x_center) + epsilon/2; % セルの中心点のx座標を計算
%y_center = x2_grid(y_center) + epsilon/2; 

% ↑で計算した平均を元に、それぞれのセルごとに矢印を書くという操作で可視化する
for i = 1:num_x1_cells
    for j = 1:num_x2_cells
        if numel(average_vecs{i,j}) > 0
            %quiver(x_center((i-1)*num_x2_cells + j), y_center((i-1)*num_x2_cells + j), average_vecs{i,j}(1)/3, average_vecs{i,j}(2)/3, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 3)
            quiver(average_vec_start_points{i, j}(1), average_vec_start_points{i, j}(2), average_vecs{i,j}(1)/3, average_vecs{i,j}(2)/3, 'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 3)
        end
    end
end
xlabel('x1');
ylabel('x2');
xlim([-0.8 0.8])
ylim([-0.7 0.7])
for i = 1:num_x1_cells
    for j = 1:num_x2_cells
        rectangle('Position', [x1_grid(i), x2_grid(j), epsilon, epsilon], 'EdgeColor', 'k'); % グリッドをプロット
    end
end
grid on;
hold off;

% セルごとの実際の平均変化ベクトルと、その変化ベクトルの始点におけるUの勾配(=F)のコサイン類似度をセルごとに算出
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

disp("平均コサイン類似度は"+cos_sim_sum / cnt + "です")


