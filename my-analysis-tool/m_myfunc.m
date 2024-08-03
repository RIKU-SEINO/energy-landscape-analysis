%% initial settings
clear all

%% params settings
syms x1 x2

% ポテンシャル関数の定義
U_sym = 200*(0.2*x1^4 + 0.4*x2^4 - 0.1*x1^2 - 0.1*x2^2);
% シミュレーションのStep幅
t_interval = 0.003;
% 揺らぎの大きさ
sigma = 2.0;
% 病気状態の平衡点
x_disease = [1/2, 1/sqrt(2)];
% 投薬後の位置
x_afterDosing = [-0.5; 0.01];

% Topic1
% シミュレーションのStep数
nPeriods_onestep = 1;
% onestepのシミュレーションの回数
total_cnt_onestep = 1;

% Topic2
% シミュレーションのStep数（Topic2）
nPeriods_manysteps = 2e5;
% シミュレーションを開始する位置
x_startPos = [0; 0];
% 分割セルの幅
gridded_interval = 0.03;

%% Topic1: Estimation probability of existence in pre-disease basin
% x_disease周りの定常状態の確率分布
disp("(START)Topic1: 病気状態平衡点周りの定常状態の確率密度関数の解析解の算出")
[pdf_average_ss, pdf_variance_ss, pdf_func_ss] = p_myfunc_linearizedPDF(x_disease, inf, sigma, U_sym);
% x_afterDosing周りに1Step時刻を進めた時の位置の確率密度関数の解析解
disp("(START)Topic1: 投薬後の位置に関する確率密度関数の解析解の算出")
[pdf_average, pdf_variance, pdf_func] = p_myfunc_linearizedPDF(x_afterDosing, t_interval, sigma, U_sym);
disp("(FINISH)Topic1: 確率密度関数の解析解の算出")

% 投薬後の位置をx_afterDosingと固定して、1step後にbasinに属している確率（解析的な確率密度関数をbasin内でintegral2で積分）
disp("(START)Topic1: basinに属する確率の計算")
stay_probability_analytic_integration = p_myfunc_stayProbabilityAnalyticIntegration(pdf_func, -Inf, 0, 0, Inf);
disp("(FINISH)Topic1: basinに属する確率の計算")


% 1step後にbasinに属している確率（実際のシミュレーション）
disp("(START)Topic1: basinに属する確率の算出")
stay_cnt = 0;
for i=1:total_cnt_onestep
    timeseries_simulation_onestep = p_myfunc_sdeSimulation(x_afterDosing, nPeriods_onestep, t_interval, sigma, U_sym);
    x1_next = timeseries_simulation_onestep(2,1);
    x2_next = timeseries_simulation_onestep(2,2);
    if (x1_next < 0) && (x2_next > 0)
        stay_cnt = stay_cnt + 1;
    end
end
stay_probability_simulation = stay_cnt / total_cnt_onestep
disp("(FINISH)Topic1: basinに属する確率の算出")


%% Topic2: Estimation of function of energy landscape
disp("(START)Topic2: ランドスケープ上のシミュレーション")
timeseries_simulation_manysteps = p_myfunc_sdeSimulation(x_startPos, nPeriods_manysteps, t_interval, sigma, U_sym);
disp("(FINISH)Topic2: ランドスケープ上のシミュレーション")

% 時系列データをdatファイルに出力
disp("(START)Topic2: 時系列データの出力")
p_myfunc_toDat(timeseries_simulation_manysteps, "samplePath.dat")
disp("(FINISH)Topic2: 時系列データの出力")

% 時系列データから遷移ベクトルデータを取得
vecs = p_myfunc_transitionVecs(timeseries_simulation_manysteps);
% 遷移ベクトルデータを各セルごとに収集
[cell_vecs, cell_vec_start_points, cell_centers] = p_myfunc_collectTransitionVecForEachCell(timeseries_simulation_manysteps, vecs, gridded_interval);
% 各セルの遷移ベクトルに関する統計情報を取得
[average_vecs, average_vec_start_points, variance_vecs, average_vec_lengths, counts] = p_myfunc_statsForEachCell(cell_vecs, cell_vec_start_points);

% コサイン類似度（uとvとwを用意し、uとv, uとwのcosine similarityを計算）
[cosine_similarity1, cosine_similarity2] = p_myfunc_cosSim(average_vecs, average_vec_start_points, cell_vec_start_points, U_sym)

% 分類
basins = p_myfunc_estimateBasin(average_vecs, cell_centers);


% グリッド
x1 = timeseries_simulation_manysteps(:,1);
x2 = timeseries_simulation_manysteps(:,2);
x1_ss_min = min(x1);
x2_ss_min = min(x2);
x1_ss_max = max(x1);
x2_ss_max = max(x2);
x1_grid = x1_ss_min:gridded_interval:x1_ss_max;
x2_grid = x2_ss_min:gridded_interval:x2_ss_max;

% 描画
disp("(START)Topic2: 描画")
%p_myfunc_drawFigure('trajectory', timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers);
%p_myfunc_drawFigure('transition_vec', timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers);
p_myfunc_drawFigure("basins", timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers)
[probability_values, energy_values] = p_myfunc_computeEnergyFromProbabilityDist(cell_vecs, sigma, x1_grid, x2_grid);
energy_values_real = p_myfunc_constructLandscape(U_sym, size(energy_values, 2), size(energy_values, 1), energy_values,  x1_grid, x2_grid);
disp("(FINISH)Topic2: 描画")