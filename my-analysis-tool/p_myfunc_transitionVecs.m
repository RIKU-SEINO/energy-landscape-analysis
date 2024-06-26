function vecs = p_myfunc_transitionVecs(timeseries)

    x1 = timeseries(:,1);
    x2 = timeseries(:,2);

    len = length(x1);
    vecs = [];
    for i = 1:len-1
        x_now = [x1(i);x2(i)];
        x_next = [x1(i+1);x2(i+1)];
        dx = x_next - x_now;
        vecs = [vecs, dx];
    end

    