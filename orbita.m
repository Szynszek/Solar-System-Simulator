%% Initialization
clc
close all
clear

%% Physical parameters
G = 6.674e-11; % [N m^2/kg^2] Gravitational Constant

bodies = {'Sun', 'Earth', 'Mars'};
body_masses = [1.989e30, 5.972e24, 6.39e23]; % [kg] Body masses

%% Initial state
% 2026-Mar-14 00:00:00.0000 TDB 
r0_sun = [-3.824460420038474e8, -8.221795612587985e8, 1.820743242961960e7]; % [m] Initial Sun position 
v0_sun = [1.212410605829541e1, 1.524019073237925e0, -2.429376929195760e-1]; % [m/s] Initial Sun velocity
y0_sun = [r0_sun, v0_sun]; % Initial Sun orbital state vector

r0_earth = [-1.480030443423478e11, 1.717356117508375e10, 1.827117222932540e7]; % [m] Initial Earth position
v0_earth = [-4.086588222194061e3, -2.968485564088994e4, 7.705410358660458e-1]; % [m/s] Initial Earth velocity
y0_earth = [r0_earth, v0_earth]; % Initial Earth orbital state vector

r0_mars = [1.754789161359977e11, -1.094239276088385e11, -6.569946055261977e9]; % [m] Initial Mars position
v0_mars = [1.366713881087808e4, 2.269063003526492e4, 1.404018870871582e2];% [m/s] Initial Mars velocity 
y0_mars = [r0_mars, v0_mars]; % Initial Mars orbital state vector

y0 = [y0_sun, y0_earth, y0_mars];

%% Simulation configuration
t_end = 31557600*2; % [s] Simulation duration
t_span = [0, t_end]; % [s] Simulation time span

%% Numerical calculations
physics_handle = @(t,y) calc_derivative(t, y, G, body_masses);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t_out, y_out] = ode45(physics_handle, t_span, y0, options);

%% Visualization
X = y_out(:, 1:6:end);
Y = y_out(:, 2:6:end);
Z = y_out(:, 3:6:end);

figure
hold on
plot3(X, Y, Z, 'LineWidth', 1.5)
legend(bodies)
axis equal
grid on
xlabel('Position X [m]')
ylabel('Position Y [m]')
zlabel('Position Z [m]')
title('Orbital Trajectory')
hold off


%% Local functions
function dydt = calc_derivative(t, y, G, body_masses)
    N = length(body_masses);
    y_formated = reshape(y, 6, N);
    r = y_formated(1:3,:);
    v = y_formated(4:6,:);
    r_i = reshape(r, 3, N, 1);
    r_j = reshape(r, 3, 1, N);
    dr = r_j-r_i;
    dist = vecnorm(dr,2,1);
    dist(dist == 0) = inf;
    m_3d = reshape(body_masses, 1, 1, N);
    da = G .* m_3d .*dr ./dist.^3;
    a = sum(da,3);

    dydt = [v; a];
    dydt = dydt(:);
end