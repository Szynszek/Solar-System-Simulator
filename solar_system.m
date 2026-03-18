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
    y0((i-1)*6+1 : i*6 ) = ephemeris_data.(body_name).state.'; 
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

y_out = integrate_PEFRL(y0, dt, num_steps, body_mu);

% %% Visualization
X = y_out(:, 1:6:end);
Y = y_out(:, 2:6:end);
Z = y_out(:, 3:6:end);

sun_idx = find(strcmp(bodies, 'Sun'), 1);
earth_idx = find(strcmp(bodies, 'Earth'), 1);
mars_idx = find(strcmp(bodies, 'Mars'), 1);

if isempty(sun_idx) || isempty(earth_idx) || isempty(mars_idx)
    error('Fatal error: Missing bodies in body vector. Check ephemeris.json file.');
end

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
plot3(X(:,sun_idx), Y(:,sun_idx), Z(:,sun_idx), '-y', 'LineWidth', 1.5, 'DisplayName', 'Sun')
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
target_idxs = [earth_idx, mars_idx];

[a_out, e_out] = calc_kepler_elements(y_out, body_mu, sun_idx, target_idxs);

figure
subplot(2,1,1)
hold on
plot(t_out/31557600, e_out(:, 1), 'b', 'LineWidth', 1.5, 'DisplayName', 'Earth')
plot(t_out/31557600, e_out(:, 2), 'r', 'LineWidth', 1.5, 'DisplayName', 'Mars')

xlabel('Time [yr]')
ylabel('Eccentricity [-]')
title('Eccentricity in years')
legend
hold off

subplot(2,1,2)
hold on
plot(t_out/31557600, a_out(:, 1)/AU, 'b', 'LineWidth', 1.5, 'DisplayName','Earth')
plot(t_out/31557600, a_out(:, 2)/AU, 'r', 'LineWidth', 1.5, 'DisplayName','Mars')

xlabel('Time [yr]')
ylabel('Semi-major axis [AU]')
title('Semi-major axis in years')
legend
hold off

% Relative and absolute error visualization
[dr_rel, dr_abs, dv_rel, dv_abs] = compare_ephemeris(y_out, bodies);
names = reshape(bodies, 1, []);

figure;

subplot(2,2,1)
plot_log_stem(dr_rel, names, 'Relative error [%]', 'Relative error of position');

subplot(2,2,2)
plot_log_stem(dv_rel, names, 'Relative error [%]', 'Relative error of velocity');

subplot(2,2,3)
plot_log_stem(dr_abs, names, 'Absolute error [m]', 'Absolute error of position');

subplot(2,2,4)
plot_log_stem(dv_abs, names, 'Absolute error [m/s]', 'Absolute error of velocity');

sgtitle(sprintf('Error analysis of PEFRL integrator | dt = %d s | T_{sim} = %.2f years', dt, t_end/31557600), 'FontWeight', 'bold');

%% Local functions

function y_out = integrate_PEFRL(y_in, dt, num_steps, body_mu)
    
    y_out = zeros(num_steps+1, size(y_in, 2));
    y_out(1,:) = y_in;
    rv = reshape(y_in, 6, []);

    r = rv(1:3,:);
    v = rv(4:6,:);

    ksi = 0.1786178958448091;
    lam = -0.2123418310626054;
    chi = -0.06626458266981849;
    
    c1 = (1-2*lam)/2 * dt;
    c2 = (1-2*(chi+ksi)) * dt;
    d1 = ksi * dt;
    d2 = chi * dt;
    d3 = lam * dt;
    
    for i = 1 : num_steps
        r1 = r + v .* d1;
        v1 = v + calc_acceleration(r1, body_mu) .* c1;
        r2 = r1 + v1 .* d2;
        v2 = v1 + calc_acceleration(r2, body_mu) .* d3;
        r3 = r2 + v2 .* c2;
        v3 = v2 + calc_acceleration(r3 , body_mu) .* d3;
        r4 = r3 + v3 .* d2;
    
        v = v3 + calc_acceleration(r4, body_mu) .* c1;
        r = r4 + v .* d1;
    
        y_out(i+1, :) = reshape([r;v], 1, []); 
    end

end

function [E_tot, dE_rel] = calc_system_energy(y_out, body_mu)
    T = size(y_out,1);
    N = size(body_mu, 3);
    G = 6.67430e-11; % [N m^2/kg^2] Gravitational Constant (CODATA)

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

function [a_out, e_out] = calc_kepler_elements(y_out, body_mu, sun_idx, target_idxs)

    T = size(y_out, 1);
    N = size(body_mu, 3);
    num_targets = length(target_idxs);

    y_3d = reshape(y_out, T, 6, N);
    y_sun = y_3d(:, :, sun_idx);
    mu_sun = body_mu(1, 1, sun_idx);

    r_rel = zeros(T, 3, num_targets);
    v_rel = zeros(T, 3, num_targets);
    mu_targets = zeros(1, 1, num_targets);

    for i = 1:num_targets
        idx = target_idxs(i);
        
        y_target_helio = y_3d(:, :, idx) - y_sun;
        r_rel(:, :, i) = y_target_helio(:, 1:3);
        v_rel(:, :, i) = y_target_helio(:, 4:6);
        
        mu_targets(1, 1, i) = mu_sun + body_mu(1, 1, idx);
    end

    r_norm = vecnorm(r_rel, 2, 2);
    v_norm = vecnorm(v_rel, 2, 2);

    a_out = squeeze((2 ./ r_norm - v_norm.^2 ./ mu_targets).^(-1));

    h_vec = cross(r_rel, v_rel, 2);
    term1 = cross(v_rel, h_vec, 2) ./ mu_targets;
    e_vec = term1 - (r_rel ./ r_norm);

    e_out = squeeze(vecnorm(e_vec, 2, 2));
end


function a = calc_acceleration(r, body_mu)

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
    a = a_newton;

end


function [dr_rel, dr_abs, dv_rel, dv_abs] = compare_ephemeris(y_in, bodies)

    raw_json = fileread('test_ephemeris.json');
    ephemeris_data = jsondecode(raw_json);
    
    N = length(bodies);
    y_test = zeros(1, 6*N);

    for i = 1:N
        body_name = bodies{i};
        
        % Verify if the simulation body exists in the test file
        if ~isfield(ephemeris_data, body_name)
            error('Fatal error: Body ''%s'' is missing in test_ephemeris.json.', body_name);
        end

        y_test((i-1)*6+1 : i*6 ) = ephemeris_data.(body_name).state.'; 
    end
    
    y_test = reshape(y_test, 6, []);
    y_in = y_in(end, :);
    y_in = reshape(y_in, 6, []);
    
    r_in = y_in(1:3, :);
    v_in = y_in(4:6, :);
    
    r_test = y_test(1:3, :);
    v_test = y_test(4:6, :);

    dr_abs = vecnorm(r_in - r_test, 2, 1); % Absolute pos error
    dv_abs = vecnorm(v_in - v_test, 2, 1); % Absolute vel error

    dr_rel = dr_abs ./ vecnorm(r_test, 2, 1) * 100; % Relative pos error
    dv_rel = dv_abs ./ vecnorm(v_test, 2, 1) * 100; % Relative vel error
end

function h = plot_log_stem(data, names, y_label, plot_title)

    data(data <= 0) = eps;
    h = stem(data, 'filled', 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', '-', 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');
    hold on; 

    xticks(1:length(names));
    xticklabels(names);
    ylabel(y_label);
    title(plot_title);
    grid on;
    set(gca, 'YScale', 'log'); 
    
    min_real = min(data);
    max_real = max(data);
    lower_limit = 10^(floor(log10(min_real)) - 0.1);
    upper_limit = 10^(ceil(log10(max_real)) + 0.1);
    ylim([lower_limit, upper_limit]);
    h.BaseValue = lower_limit; 
    hold off;
end