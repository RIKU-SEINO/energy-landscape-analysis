function energy_values_real = p_myfunc_realEnergyValues(U_sym, interval_x1, interval_x2)
    % シンボリック変数の定義
    syms x1 x2

    % 数値計算用の関数に変換
    U = matlabFunction(U_sym, 'Vars', [x1, x2]);
  
    % x1とx2の範囲を設定
    x1_range = linspace(-0.8, 0.8, interval_x1);
    x2_range = linspace(-0.7, 0.7, interval_x2);
  
    % メッシュグリッドを作成
    [X1, X2] = meshgrid(x1_range, x2_range);
  
    % ポテンシャル関数の値を計算
    U_vals = U(X1, X2);
  
    % セル配列を作成し、エネルギー値を格納
    energy_values_real = cell(size(U_vals));
    for i = 1:size(U_vals, 1)
        for j = 1:size(U_vals, 2)
            energy_values_real{i, j} = U_vals(i, j);
        end
    end