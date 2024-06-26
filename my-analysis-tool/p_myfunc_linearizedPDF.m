function result = p_myfunc_linearizedPDF(x0, t_interval, sigma, U_sym)

    % シンボリック変数の定義
    syms x1 x2

    x = [x1;x2];

    % 勾配の計算（各変数に対する偏微分）
    grad_U = gradient(U_sym, [x1, x2]);
    grad_U_func = matlabFunction(grad_U, 'Vars', {x1, x2});

    % ヘッセ行列の計算
    hessian_U = hessian(U_sym, [x1, x2]);
    hessian_U_func = matlabFunction(hessian_U, 'Vars', {x1, x2});

    % OU過程のパラメータ
    mu = - inv(hessian_U_func(x0(1), x0(2))) * grad_U_func(x0(1), x0(2));
    theta = hessian_U_func(x0(1), x0(2));

    [V, theta_tilde] = eig(theta);
    lambda = diag(theta_tilde);

    diag_elements = (1 ./ (2 .* lambda)) .* (1 - exp(-2 .* lambda .* t_interval));
    D_new = diag(diag_elements);

    average = x0 + (eye(length(x0)) - expm(-theta*t_interval))*mu;
    variance = sigma.^2 * V * D_new * V.';

    result = containers.Map();
    result('average') = average;
    result('variance') = variance;
    result('pdf') = (1 / (2 * pi * sqrt(det(variance)))) * expm(-0.5 * (x - average).' * inv(variance) * (x - average));