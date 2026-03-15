%% Initialization
clc
close all
clear

%% Physical parameters
G = 6.674e-11; % [N m^2/kg^2] Gravitational Constant
AU = 1.495978707e11; % [m] Astronomical unit 


bodies = {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
body_masses = [1.989e30, 3.301e23, 4.867e24, 6.045e24, 6.39e23, 1.898e27, 5.683e26, 8.680e25, 1.024e26]; % [kg] Body masses (Earth+Moon mass)
body_masses = reshape(body_masses, 1, 1, []);

%% Initial state
% 2026-Mar-14 00:00:00.0000 TDB 
r0_sun = [-382446042.0038474, -822179561.2587985, 18207432.4296196];
v0_sun = [12.124106058295409, 1.524019073237925, -0.242937692919576];
y0_sun = [r0_sun, v0_sun];

r0_mercury = [-59484010545.39254, -14138364260.07398, 4350629859.97398];
v0_mercury = [602.0470130047652, -45426.86225273868, -3766.897484204248];
y0_mercury = [r0_mercury, v0_mercury];

r0_venus = [92153321309.08293, 55328018669.09437, -4549636525.406504];
v0_venus = [-18265.144344563003, 29790.35226168731, 1463.605242061472];
y0_venus = [r0_venus, v0_venus];

% Changed to Earth-Moon barycenter
r0_earth = [-148001259534.1553, 17169091157.38058, 17947884.04517062];
v0_earth = [-4075.6727790927607, -29679.940956888782, 1.560242767398634];
y0_earth = [r0_earth, v0_earth];

r0_mars = [175478916135.9977, -109423927608.8385, -6569946055.261977];
v0_mars = [13667.138810878081, 22690.63003526492, 140.4018870871582];
y0_mars = [r0_mars, v0_mars];

r0_jupiter = [-330177474985.5681, 709957964651.7996, 4444301040.555298];
v0_jupiter = [-12003.139051740809, -4891.552424424635, 288.88157796239625];
y0_jupiter = [r0_jupiter, v0_jupiter];

r0_saturn = [1415642229227.292, 97616615972.22827, -58062166148.78422];
v0_saturn = [-1196.585897214322, 9616.828736369385, -120.1764782487014];
y0_saturn = [r0_saturn, v0_saturn];

r0_uranus = [1440633774263.0889, 2531658218329.792, -9261249460.201502];
v0_uranus = [-5968.7510279623775, 3050.593482898619, 88.60253443041711];
y0_uranus = [r0_uranus, v0_uranus];

r0_neptune = [4467414306173.692, 110809131133.7783, -105238080957.0405];
v0_neptune = [-171.2848260947789, 5465.766243943168, -108.0812461311418];
y0_neptune = [r0_neptune, v0_neptune];

y0 = [y0_sun, y0_mercury, y0_venus, y0_earth, y0_mars, y0_jupiter, y0_saturn, y0_uranus, y0_neptune];


%% Simulation configuration
t_end = 31557600*30; % [s] Simulation duration
t_span = [0, t_end]; % [s] Simulation time span

%% Numerical calculations
dt = 3600;
num_steps = round(t_end / dt);
t_out = linspace(0, t_end, num_steps+1);

N = length(bodies);
y_out = zeros(num_steps+1, 6*N);
y_out(1,:) = y0;
y_formatted = reshape(y0, 6, N);

r = y_formatted(1:3,:);
v = y_formatted(4:6,:);
a = calc_acceleration(r, v, G, body_masses);

for i = 1 : num_steps
    v_half = v + a .* (dt/2);
    r = r + v_half .* dt;
    a = calc_acceleration(r, v_half, G, body_masses);
    v = v_half + a .* (dt/2);
    y_out(i+1,:) = reshape([r;v],1 , []);
end

%% Visualization
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
[E_tot, dE_rel] = calc_system_energy(y_out, G, body_masses);
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
[a_out, e_out] = calc_kepler_elements(y_out, G, body_masses);

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

function [E_tot, dE_rel] = calc_system_energy(y_out, G, body_masses)
    T = size(y_out,1);
    N = size(body_masses, 3);

    y_3d = reshape(y_out, T, 6, N); % Planet state per page

    r = y_3d(:,1:3,:); % Only positions
    v = y_3d(:,4:6,:); % Only velocity
    v2 = v.^2;
    
    E_k = 0.5 .* sum(v2,2) .* body_masses;
    E_k = sum(E_k,3);
    E_k = E_k(:);
    
    r_i = reshape(r, T, 3, N, 1); % 4d tensor
    r_j = reshape(r, T, 3, 1, N); % 4d tensor
    dr = r_j-r_i; % 4d tensor (Time, position, body_i, body_j)

    dist = vecnorm(dr, 2, 2);
    dist(dist==0) = inf;

    m_i = reshape(body_masses, 1, 1, N, 1);
    m_j = reshape(body_masses, 1, 1, 1, N);
    
    E_p = -G .* (m_i .* m_j) ./ dist;

    E_p = sum(E_p, [3, 4]) / 2;
    E_p = E_p(:);
    
    E_tot = E_k + E_p;

    dE_rel = (E_tot - E_tot(1)) ./ abs(E_tot(1)); % Relative energy error
end

function [a_out, e_out] = calc_kepler_elements(y_out, G, body_masses)

    N = size(body_masses, 3);
    T = size(y_out, 1);
    y_3d = reshape(y_out, T, 6, N);
    y_sun = y_out(:, 1:6);
    y_sun = reshape(y_sun, T, 6, 1);
    
    y_helio = y_3d - y_sun;
    r_rel = y_helio(:, 1:3, 2:N); % without page with sun
    v_rel = y_helio(:, 4:6, 2:N); % without page with sun
    
    m_sun = body_masses(:,:,1);
    m_planets = body_masses(:,:, 2:N);
    mu = G .* (m_planets + m_sun);
    
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

function a = calc_acceleration(r, v, G, body_masses)

    c = 299792458; % [m/s] Speed of light
    N = size(body_masses, 3);
    
    r_i = reshape(r, 3, N, 1);
    r_j = reshape(r, 3, 1, N);
    dr = r_j-r_i;
    dist = vecnorm(dr,2,1);
    dist(dist == 0) = inf;
    da = G .* body_masses .*dr ./dist.^3;
    a_newton = sum(da,3);
    a_newton = squeeze(a_newton);
    
    m_sun = body_masses(:,:,1);
    r_sun = r(:, 1);
    v_sun = v(:, 1);

    r_rel = r - r_sun;
    v_rel = v - v_sun;
    
    r_norm = vecnorm(r_rel, 2, 1);
    v_norm = vecnorm(v_rel, 2, 1);

    r_norm(1) = inf;

    
    mu = G * m_sun;
    r_dot_v = sum((r_rel .* v_rel), 1);
    

    a_GR = mu ./ (c^2 .* r_norm.^3) .* ((4 * mu ./ r_norm - v_norm.^2) .* r_rel + 4 .* (r_dot_v) .* v_rel); % Schwarzschild solution

    a = a_newton + a_GR;

end