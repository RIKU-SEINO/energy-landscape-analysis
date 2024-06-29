function result = p_myfunc_sdeSimulation(x0, nPeriods, t_interval, sigma, U_sym)

    % シンボリック変数の定義
    syms x1 x2
    
    % 勾配の計算（各変数に対する偏微分）
    grad_U = gradient(U_sym, [x1, x2]);
    
    % 勾配の数値関数の定義
    grad_U_func = matlabFunction(grad_U, 'Vars', {x1, x2});
    
    % ドリフト項
    drift = @(t,x)-grad_U_func(x(1), x(2));
    
    % 拡散項
    diffusion = @(t,x) sigma*eye(2);
    
    % sdeオブジェクトの作成
    obj = sde(drift, diffusion, 'StartState', x0);
    
    % Simulate the system
    [result,~] = simByEuler(obj, nPeriods, 'DeltaTime', t_interval);

    % 打ち切りする条件を設定（ここの値を超えたらそれ以降は発散するので非現実的）
    threshold = 5;

    % 行列を走査して条件を満たしたら打ち切る
    for i = 1:size(result, 1)
        if any(abs(result(i, :)) > threshold)
            disp("打ち切りを行いました")
            result = result(1:i, :);
            break;
        end
    end
end