%% initial settings
clear all


%% params settings
syms x1 x2

% ポテンシャル関数の定義
U_sym = 200*(0.2*x1^4 + 0.4*x2^4 - 0.1*x1^2 - 0.1*x2^2);
% シミュレーションのStep幅
t_interval = 0.01;
% 揺らぎの大きさ
sigma = 1.6;

% Topic1
% シミュレーションのStep数
nPeriods_onestep = 1;
% 投薬後の位置
x_afterDosing = [-0.5; 0.01];
% onestepのシミュレーションの回数
total_cnt_onestep = 1;

% Topic2
% シミュレーションのStep数（Topic2）
nPeriods_manysteps = 5e4;
% シミュレーションを開始する位置
x_startPos = [0; 0];
% 分割セルの幅
gridded_interval = 0.05;


%% Topic1: Estimation probability of existence in pre-disease basin
% x_afterDosing周りに1Step時刻を進めた時の位置の確率密度関数の解析解
disp("(START)Topic1: 確率密度関数の解析解の算出")
[pdf_average, pdf_variance, pdf_func] = p_myfunc_linearizedPDF(x_afterDosing, t_interval, sigma, U_sym);
disp("(FINISH)Topic1: 確率密度関数の解析解の算出")


% 1step後にbasinに属している確率（解析的な確率密度関数をbasin内でintegral2で積分）
disp("(START)Topic1: basinに属する確率の計算")
stay_probability_analytic_integration = p_myfunc_stayProbabilityAnalyticIntegration(pdf_func, -Inf, 0, 0, Inf)
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

% コサイン類似度（uとw）
cosine_similarity = p_myfunc_cosSim(average_vecs, average_vec_start_points, U_sym)

% コサイン類似度（uとv）

% 分類
basins = p_myfunc_estimateBasin(average_vecs, cell_centers);

% 時系列データのパスを描画
disp("(START)Topic2: 描画")
%p_myfunc_drawFigure('trajectory', timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers);
p_myfunc_drawFigure('transition_vec', timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers);
p_myfunc_drawFigure("basins", timeseries_simulation_manysteps, average_vecs, average_vec_start_points, U_sym, gridded_interval, basins, cell_centers)
disp("(FINISH)Topic2: 描画")




