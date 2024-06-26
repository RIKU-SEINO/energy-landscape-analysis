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

% Simulate the system
epsilon = 0.05;% granularity of grid
nPeriods = 1;% simulation step length
dt = 0.01; % time step width of simulation
stationaryStartIndex = 1; % index to start to record
nSimulationByCell = 3;
sigma = 0.8;


% U(x1,x2) = 200*(0.2*x1.^4 + 0.4*x2.^4 -0.1*x1.^2 -0.1*x2.^2);
% gradient of U = [(4*x1.^3)/5 - x1/5, (8*x2^3)/5 - x2/5]

% Define the system of equations
F = @(t,X) 200*[ - (4*X(1).^3)/5 + X(1)/5; -(8*X(2).^3)/5 + X(2)/5 ];
G = @(t,X) [sigma, 0; 0, sigma];

% grid params
x1_ss_min = -1;
x2_ss_min = -1;
x1_ss_max = 1;
x2_ss_max = 1;
x1_grid = x1_ss_min:epsilon:x1_ss_max;
x2_grid = x2_ss_min:epsilon:x2_ss_max;
num_x1_cells = length(x1_grid);
num_x2_cells = length(x2_grid);
cell_vecs = cell(num_x1_cells, num_x2_cells);%各セル内に存在するinitialStateからの変化ベクトルをセルごとに全てを格納、そのための初期化
cell_vecs_gradient = cell(num_x1_cells, num_x2_cells);%各セル内に存在するinitialStateにおけるgradientベクトルをセルごとに全てを格納、そのための初期化
cell_vec_start_points = cell(num_x1_cells, num_x2_cells);% 各セル内に存在する変化ベクトルの始点をセルごとに全てを格納、そのための初期化
% num_x1_cells × num_x2_cellsで区切られたgridによる各々のセルについてnSimulationByCell回プロット&シミュレーションを行い、そのデータを取得する。
for i = 1:num_x1_cells-1
    for j = 1:num_x2_cells-1
        x1_ss_cell = x1_grid(:,i:i+1);
        x2_ss_cell = x2_grid(:,i:i+1);
        x1_ss_min_cell = x1_ss_cell(1);
        x1_ss_max_cell = x1_ss_cell(2);
        x2_ss_min_cell = x2_ss_cell(1);
        x2_ss_max_cell = x2_ss_cell(2);
        x1_samples = x1_ss_min_cell + (x1_ss_max_cell - x1_ss_min_cell) * rand(nSimulationByCell, 1);
        x2_samples = x2_ss_min_cell + (x2_ss_max_cell - x2_ss_min_cell) * rand(nSimulationByCell, 1);
        for sample_index = 1:nSimulationByCell
            initialState = [x1_samples(sample_index);x2_samples(sample_index)];
            obj = sde(F, G,'StartState', initialState);
            [S,T] = simByEuler(obj, nPeriods, 'DeltaTime', dt);
            x1 = S(:,1);
            x2 = S(:,2);
            x1_stationary = x1(stationaryStartIndex:end);
            x2_stationary = x2(stationaryStartIndex:end);
            vec = [x1_stationary(2)-x1_stationary(1);x2_stationary(2)-x2_stationary(1)];
            if isempty(cell_vecs{i, j})
                cell_vecs{i, j} = vec;
                cell_vec_start_points{i, j} = initialState;
                cell_vecs_gradient{i ,j} = F(0, initialState);
            else
                cell_vecs{i, j} = [cell_vecs(i, j), vec];
                cell_vec_start_points{i, j} = [cell_vec_start_points{i, j}, initialState];
                cell_vecs_gradient{i ,j} = [cell_vecs_gradient{i, j}, F(0, initialState)];
            end
        end
    end
end



