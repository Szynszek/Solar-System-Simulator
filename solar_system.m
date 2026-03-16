%% Initialization
clc
close all
clear

raw_json = fileread('ephemeris.json');
ephemeris_data = jsondecode(raw_json);
bodies = fieldnames(ephemeris_data);
N = length(bodies);
body_mu = zeros(1, 1, N);
y0 = zeros(1, 6*N);

for i = 1:N
    body_name = bodies{i};
    body_mu(1, 1, i) = ephemeris_data.(body_name).mu;
    y0((i-1)*6+1 : i*6 ) = ephemeris_data.(body_name).state'; 
end


%% Physical parameters
AU = 1.495978707e11; % [m] Astronomical unit 

%% Simulation configuration
t_end = 31557600*30; % [s] Simulation duration
t_span = [0, t_end]; % [s] Simulation time span

%% Numerical calculations
dt = 1800; % [s] Simulation time step size
num_steps = round(t_end / dt);
t_out = linspace(0, t_end, num_steps+1);

y_out = zeros(num_steps+1, 6*N);
y_out(1,:) = y0;
y_formatted = reshape(y0, 6, N);

r = y_formatted(1:3,:);
v = y_formatted(4:6,:);
a = calc_acceleration(r, v, body_mu);

for i = 1 : num_steps
    v_half = v + a .* (dt/2);
    r = r + v_half .* dt;
    a = calc_acceleration(r, v_half, body_mu);
    v = v_half + a .* (dt/2);
    y_out(i+1,:) = reshape([r;v],1 , []);
end

% %% Visualization
X = y_out(:, 1:6:end);
Y = y_out(:, 2:6:end);
Z = y_out(:, 3:6:end);

% Orbital visualization
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

% Sun trajectory visualization
figure
hold on
plot3(X(:,1), Y(:,1), Z(:,1), '-y', 'LineWidth', 1.5, 'DisplayName','Sun')
plot3(0,0,0, '*r', 'DisplayName','Barycenter')
legend
axis equal
grid on
xlabel('Position X [m]')
ylabel('Position Y [m]')
zlabel('Position Z [m]')
title('Sun orbital trajectory')
hold off

% Total energy and error visualization
[E_tot, dE_rel] = calc_system_energy(y_out, body_mu);
figure
plot(t_out, E_tot)
xlabel('Time [s]')
ylabel('Total energy [J]')
title('Total energy in time')
grid on

figure
plot(t_out, dE_rel)
xlabel('Time [s]')
ylabel('Relative energy error [-]')
title('Relative energy error in time')
grid on

% Eccentricity and Semi-major axis change visualization
[a_out, e_out] = calc_kepler_elements(y_out, body_mu);

subplot(2,1,1)
hold on
plot(t_out/31557600, e_out(:, 3), 'b', 'LineWidth',1.5, 'DisplayName','Earth')
plot(t_out/31557600, e_out(:, 4), 'r', 'LineWidth',1.5, 'DisplayName','Mars')

xlabel('Time [yr]')
ylabel('Eccentricity [-]')
title('Eccentricity in years')
legend
hold off

subplot(2,1,2)
hold on
plot(t_out/31557600, a_out(:, 3)/AU, 'b', 'LineWidth',1.5, 'DisplayName','Earth')
plot(t_out/31557600, a_out(:, 4)/AU, 'r', 'LineWidth',1.5, 'DisplayName','Mars')

xlabel('Time [yr]')
ylabel('Semi-major axis [AU]')
title('Semi-major axis in years')
legend
hold off

%% Local functions

function [E_tot, dE_rel] = calc_system_energy(y_out, body_mu)
    T = size(y_out,1);
    N = size(body_mu, 3);
    G = 6.674e-11; % [N m^2/kg^2] Gravitational Constant

    y_3d = reshape(y_out, T, 6, N); % Planet state per page

    r = y_3d(:,1:3,:); % Only positions
    v = y_3d(:,4:6,:); % Only velocity
    v2 = v.^2;
    
    E_k = 0.5 .* sum(v2,2) .* (body_mu ./ G);
    E_k = sum(E_k,3);
    E_k = E_k(:);
    
    r_i = reshape(r, T, 3, N, 1); % 4d tensor
    r_j = reshape(r, T, 3, 1, N); % 4d tensor
    dr = r_j-r_i; % 4d tensor (Time, position, body_i, body_j)

    dist = vecnorm(dr, 2, 2);
    dist(dist==0) = inf;

    mu_i = reshape(body_mu, 1, 1, N, 1);
    mu_j = reshape(body_mu, 1, 1, 1, N);
    
    E_p = -(mu_i .* mu_j) ./ (dist .* G);

    E_p = sum(E_p, [3, 4]) / 2;
    E_p = E_p(:);
    
    E_tot = E_k + E_p;

    dE_rel = (E_tot - E_tot(1)) ./ abs(E_tot(1)); % Relative energy error
end

function [a_out, e_out] = calc_kepler_elements(y_out, body_mu)

    N = size(body_mu, 3);
    T = size(y_out, 1);
    y_3d = reshape(y_out, T, 6, N);
    y_sun = y_out(:, 1:6);
    y_sun = reshape(y_sun, T, 6, 1);
    
    y_helio = y_3d - y_sun;
    r_rel = y_helio(:, 1:3, 2:N); % without page with sun
    v_rel = y_helio(:, 4:6, 2:N); % without page with sun
    
    mu = body_mu(:,:,2:N) + body_mu(:,:,1); % mu = G * (m_p + M_s)
    
    r_norm = vecnorm(r_rel, 2, 2);
    v_norm = vecnorm(v_rel, 2, 2);
    
    a_out = (2./r_norm - v_norm.^2 ./ mu).^(-1);
    a_out = squeeze(a_out);
    
    h_vec = cross(r_rel, v_rel, 2);
    
    term1 = cross(v_rel, h_vec, 2) ./ mu;
    
    e_vec = term1 - (r_rel./r_norm);
    
    e_out = vecnorm(e_vec, 2, 2);
    e_out = squeeze(e_out);
    
end

function a = calc_acceleration(r, v, body_mu)

    c = 299792458; % [m/s] Speed of light
    N = size(body_mu, 3);
    
    % Newton term
    r_i = reshape(r, 3, N, 1);
    r_j = reshape(r, 3, 1, N);
    dr = r_j-r_i;
    dist = vecnorm(dr,2,1);
    dist(dist == 0) = inf;
    da = body_mu .*dr ./dist.^3;
    a_newton = sum(da,3);
    a_newton = squeeze(a_newton);
    
    % Schwarzschild solution term
    mu = body_mu(:,:,1);
    r_sun = r(:, 1);
    v_sun = v(:, 1);

    r_rel = r - r_sun; % 3 x N-1 
    v_rel = v - v_sun;
    
    r_norm = vecnorm(r_rel, 2, 1);
    v_norm = vecnorm(v_rel, 2, 1);

    r_norm(1) = inf;
    
    r_dot_v = sum((r_rel .* v_rel), 1);
    
    a_GR = mu ./ (c^2 .* r_norm.^3) .* ((4 * mu ./ r_norm - v_norm.^2) .* r_rel + 4 .* (r_dot_v) .* v_rel); % Schwarzschild solution

    a = a_newton + a_GR;

end