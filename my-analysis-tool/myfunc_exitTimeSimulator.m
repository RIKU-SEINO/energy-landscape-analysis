clear all
close all

% Parameters
num_simulations = 5000; % Number of simulations
dt = 0.01; % Time step size
initial_state = [-0.1; 1/(2*sqrt(2))]; % Initial state in well A (e.g., 2nd quadrant)

% Define potential function and its gradient
U = @(x, y) 200 * (0.2 * x.^4 + 0.4 * y.^4 - 0.1 * x.^2 - 0.1 * y.^2);
F = @(t,X) 200*[ - (4*X(1).^3)/5 + X(1)/5; -(8*X(2).^3)/5 + X(2)/5 ];
G = @(t,X) [ 0.8, 0; 0, 0.8];
sigma = 0.8;
grad2_U = @(t, X) 200 * [1/5 - 12*X(1).^2 / 5, 0; 0, 1/5 - 24*X(2).^2 / 5];

% define energy landscape function
figure(1)
x = linspace(-0.8,0.8);
y = linspace(-0.7,0.7);
[x,y] = meshgrid(x,y);
z = U(x,y);
contour(x, y, z, 50);
title('Surface Plot');
% Show the colorbar for reference
colorbar;
hold on;

% Create the sde object
obj = sde(F, G, 'StartState', initial_state);

% Initialize transition counters for each basin
transitions = [0, 0, 0, 0];

% Run simulations
for sim = 1:num_simulations
    % Simulate the system
    [S,T] = simByEuler(obj, 1, 'DeltaTime', dt);
    X = S(1,:);
    X_next = S(2,:);

    if (X_next(1)*X(1) < 0) && (X_next(2)*X(2) > 0)%transition to first quadrant
        transitions(1) = transitions(1) + 1; 
    elseif (X_next(1)*X(1) < 0) && (X_next(2)*X(2) < 0)%transition to fourth quadrant
        transitions(4) = transitions(4) + 1; 
    elseif (X_next(1)*X(1) > 0) && (X_next(2)*X(2) < 0)%transition to third quadrant
        transitions(3) = transitions(3) + 1; 
    else
        transitions(2) = transitions(2) + 1; 
    end
    %quiver(X(1), X(2), X_next(1) -  X(1), X_next(2) -  X(2), 0, 'r', 'LineWidth', 1, 'MaxHeadSize', 0.5);
    plot(X_next(1), X_next(2), 'b*')
    hold on;
end
hold on;
plot(initial_state(1), initial_state(2), 'ro');


% Calculate the transition rates for each basin
transition_rates = transitions / num_simulations;

% Display the results
fprintf('Estimated transition rates:\n');
fprintf('To first quadrant: %f\n', transition_rates(1));
fprintf('To second quadrant: %f\n', transition_rates(2));
fprintf('To third quadrant: %f\n', transition_rates(3));
fprintf('To fourth quadrant: %f\n', transition_rates(4));

% Analytical transition rate for comparison (assuming symmetry and example values)
A = [prefactor([0;1/(2*sqrt(2))], initial_state);
          0;
          prefactor([-1/2;0], initial_state);
          prefactor([0;0], initial_state)];
delta_U = [U(0,1/(2*sqrt(2)))-U(initial_state(1),initial_state(2));
           0;
           U(-1/2,0)-U(initial_state(1),initial_state(2));
           U(0,0)-U(initial_state(1),initial_state(2))]; 
analytical_transition_rate = [A(1) * exp(- 2*delta_U(1) / sigma^2);
                              A(2) * exp(- 2*delta_U(2) / sigma^2);
                              A(3) * exp(- 2*delta_U(3) / sigma^2);
                              A(4) * exp(- 2*delta_U(4) / sigma^2)];
analytical_transition_rate = analytical_transition_rate / sum(analytical_transition_rate);                   
fprintf('Analytical transition rate: %f\n', analytical_transition_rate);

function out = prefactor(saddle_point, eq_point)
    grad2_U = @(t, X) 200 * [1/5 - 12*X(1).^2 / 5, 0; 0, 1/5 - 24*X(2).^2 / 5];
    hessian_saddle = grad2_U(0, saddle_point);
    hessian_eq = grad2_U(0, eq_point);

    det_hessian_saddle = det(hessian_saddle);
    det_hessian_eq = det(hessian_eq);

    eig_values = eig(hessian_saddle);
    for i = 1:length(eig_values)
        eig_value = eig_values(i);
        if eig_value < 0
            break
        else
            eig_value = 0;
        end
    end

    out = abs((eig_value/(2*pi)) * sqrt(det_hessian_eq / det_hessian_saddle));
end
