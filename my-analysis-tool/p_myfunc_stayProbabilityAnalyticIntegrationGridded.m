function result = p_myfunc_stayProbabilityAnalyticMonteCarlo(pdf, vertices)

    eps = 1e-5;
    infty_value = 1e3;

    for i = 1:length(vertices)
        if abs(vertices(i,1)) > infty_value
            vertices(i,1) = sign(vertices(i,1))*infty_value;
        end
        if abs(vertices(i,2)) > infty_value
            vertices(i,2) = sign(vertices(i,2))*infty_value;
        end
    end
    disp(vertices)
   
    for i = 1:length(vertices)
        [edge_vec_x, edge_vec_y] = edge_vec_cal(vertices, i);
        if abs(edge_vec_x) < eps
            vertices(i,1) = vertices(i,1) + 1e-2*std(vertices(:,1))*rand(1);
        end
        if abs(edge_vec_y) < eps
            vertices(i,2) = vertices(i,2) + 1e-2*std(vertices(:,2))*rand(1);
        end
    end
       
    x = vertices(:,1);
    y = vertices(:,2);
    disp(vertices(2,1))

    pdf_func = matlabFunction(pdf);
    
    if isrow(x)
        x=x';
    end
    
    if isrow(y)
        y=y';
    end

    % delete duplicate verterces
    mypolygon = unique([x,y],'rows','stable');
    x = mypolygon(:,1);
    y = mypolygon(:,2);

    %[x,y]=poly2cw(x,y);
    [xmin,ind1] = min(x);
    x = circshift(x,-(ind1-1));
    y = circshift(y,-(ind1-1));
    [xmax,ind2] = max(x);

    if y(2) > y(end)
        up = 1:ind2;
        down = [ind2:length(x),1];
    else
        down = 1:ind2;
        up = [ind2:length(x),1];
    end

    ymin = @(xx)interp1(x(down),y(down),xx) ;%regional lower bound
    ymax = @(xx)interp1(x(up),y(up),xx) ;%regional  upper bound
    result = integral2(pdf_func,xmin,xmax,ymin,ymax);
end

function [x, y] = edge_vec_cal(vertices, idx)
    vertex = vertices(idx,:);
    if idx == length(vertices)
        vertex_next = vertices(1,:);
    else
        vertex_next = vertices(idx+1,:);
    end
    vec = vertex_next - vertex;
    x = vec(1);
    y = vec(2);
end


