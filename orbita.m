%% Inicjalizacja
clc
close all
clear

%% Parametry fizyczne
G = 6.674e-11; % [N m^2/kg^2] Stała grawitacji

planety = {'Słońce', 'Ziemia', 'Mars'};
masy_planet = [1.989e30, 5.972e24, 6.39e23]; % [kg] Masy Ciał

%% Stan początkowy
r0_slonce = [0, 0, 0]; % [m] Początkowe położenie Słońca
v0_slonce = [0, 0, 0]; % [m/s] Początkowa prędkość Słońca
y0_slonce = [r0_slonce, v0_slonce]; % Wektor stanu początkowego Słońca

r0_ziemia = [0, 1.498e11, 0]; % [m] Początkowe położenie Ziemi
v0_ziemia = [-29.78e3 * 0.70, 0, 0]; % [m/s] Początkowa prędkość Ziemi
y0_ziemia = [r0_ziemia, v0_ziemia]; % Wektor stanu początkowego Ziemi

r0_mars = [0, 2.279e11, 0]; % [m] Początkowe położenie Marsa
v0_mars = [-24.07e3, 0, 0];% [m/s] Początkowa prędkość Marsa 
y0_mars = [r0_mars, v0_mars]; % Wektor stanu początkowego Marsa

y0 = [y0_slonce, y0_ziemia, y0_mars];

%% Ustawienia symulacji
t_end = 31557600*2; % [s] Czas trwania symulacji
t_span = [0, t_end]; % [s] Przedział czasu symulacji

%% Obliczenia numeryczne
uchwyt_fizyki = @(t,y) oblicz_pochodne(t, y, G, masy_planet);
opcje = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t_out, y_out] = ode45(uchwyt_fizyki, t_span, y0, opcje);

%% Wizualizacja
X = y_out(:, 1:6:end);
Y = y_out(:, 2:6:end);
Z = y_out(:, 3:6:end);

figure
hold on
plot3(X, Y, Z, 'LineWidth', 1.5)
legend(planety)
axis equal
grid on
xlabel('Pozycja X [m]')
ylabel('Pozycja Y [m]')
zlabel('Pozycja Z [m]')
title('Trajektoria orbitalna')
hold off


%% Funkcje
function dydt = oblicz_pochodne(t, y, G, masy_planet)
    N = length(masy_planet);
    y_formated = reshape(y, 6, N);
    r = y_formated(1:3,:);
    v = y_formated(4:6,:);
    r_i = reshape(r, 3, N, 1);
    r_j = reshape(r, 3, 1, N);
    dr = r_j-r_i;
    dist = vecnorm(dr,2,1);
    dist(dist == 0) = inf;
    m_3d = reshape(masy_planet, 1, 1, N);
    da = G .* m_3d .*dr ./dist.^3;
    a = sum(da,3);

    dydt = [v; a];
    dydt = dydt(:);
end