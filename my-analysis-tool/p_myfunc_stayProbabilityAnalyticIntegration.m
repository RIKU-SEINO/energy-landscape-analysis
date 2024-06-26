function result = p_myfunc_stayProbabilityAnalyticIntegration(pdf, lowerLimitX1, upperLimitX1, lowerLimitX2, upperLimitX2)
    syms x1 x2
    pdf_func = matlabFunction(pdf, 'vars', {x1, x2});
    result = integral2(pdf_func, lowerLimitX1, upperLimitX1, lowerLimitX2, upperLimitX2); % 積分領域を四角形として計算しているが、今後多角形に一般化する必要ありそうかも
