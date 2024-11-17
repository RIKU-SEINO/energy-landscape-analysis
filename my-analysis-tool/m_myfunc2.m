% パラメータの設定（適宜調整してください）
params.mu1 = 0.1;
params.mu2 = 0.041;
params.alpha = 0.052;
params.r = 0.0032;
params.beta = 1.1e-8;
params.p1 = 1.25e-7;
params.p2 = 0.285e-7;
params.p3 = 1.1e-7;
params.p4 = 0.12e-7;
params.p5 = 0.345e-9;

% 初期条件とシミュレーションの時間範囲
Y0 = [1; 1; 1; 1];  % 初期条件 [B(0), E(0), Ti(0), Tu(0)]
tspan = [0 1e4];     % 時間範囲

% 微分方程式の数値解を求める
[t, Y] = ode45(@(t, Y) system(t, Y, params), tspan, Y0);

% 解の取り出し
B = Y(:,1);
E = Y(:,2);
Ti = Y(:,3);
Tu = Y(:,4);

% 3Dプロット (E, Ti, Tu の軌跡)
figure;
plot3(E, Ti, Tu, '-o', 'LineWidth', 1.5);
xlabel('E(t)');
ylabel('T_i(t)');
zlabel('T_u(t)');
title('3D Trajectory of E, T_i, T_u');
grid on;

% 微分方程式系の定義
function dYdt = system(t, Y, params)
    B = Y(1);
    E = Y(2);
    Ti = Y(3);
    Tu = Y(4);

    dB_dt = -params.mu1 * B - params.p1 * E * B - params.p2 * B * Tu;
    dE_dt = -params.mu2 * E + params.alpha * Ti + params.p4 * E * B - params.p5 * E * Ti;
    dTi_dt = -params.p3 * E * Ti + params.p2 * B * Tu;
    dTu_dt = -params.p2 * B * Tu + params.r * (1 - params.beta * Tu) * Tu;

    dYdt = [dB_dt; dE_dt; dTi_dt; dTu_dt];
end
