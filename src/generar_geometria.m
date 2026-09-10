%% A.1 – GENERACIÓN DE LA GEOMETRÍA Y PANELADO
%% Perfil NACA 6409 – Distribución coseno

clear; clc; close all;

%% 1. PARÁMETROS DEL PERFIL (NACA 6409)
f_max = 0.06;   % curvatura máxima (6%)
x_fmax = 0.40;   % posición de curvatura máxima
t_max = 0.09;   % espesor máximo (9%)
c = 1.0;   % cuerda
N_half = 50;   % paneles por cara (total N = 100)

%% 2. DISTRIBUCIÓN COSENO DE LOS NODOS EN x
beta = linspace(0, pi, N_half+1);
x_dist = 0.5 * (1 - cos(beta));

%% 3. LÍNEA DE CURVATURA MEDIA Y DERIVADA
zc = zeros(size(x_dist));
dzc_dx = zeros(size(x_dist));

for i = 1:length(x_dist)
    x = x_dist(i);

    if x < x_fmax
        zc(i) = (f_max/x_fmax^2)*(2*x_fmax*x - x^2);
        dzc_dx(i) = (2*f_max/x_fmax^2)*(x_fmax - x);
    else
        zc(i) = (f_max/(1-x_fmax)^2)*((1-2*x_fmax) + 2*x_fmax*x - x^2);
        dzc_dx(i) = (2*f_max/(1-x_fmax)^2)*(x_fmax - x);
    end
end

theta_geo = atan(dzc_dx);

%% 4. DISTRIBUCIÓN DE ESPESOR
ze = 5*t_max * (0.2969*sqrt(x_dist) - 0.1260*x_dist ...
    - 0.3516*x_dist.^2 + 0.2843*x_dist.^3 - 0.1015*x_dist.^4);

%% 5. COORDENADAS EXTRADÓS E INTRADÓS
xu = x_dist - ze .* sin(theta_geo);
zu = zc + ze .* cos(theta_geo);

xl = x_dist + ze .* sin(theta_geo);
zl = zc - ze .* cos(theta_geo);

%% 6. CONTORNO CERRADO DEL PERFIL
X = [flip(xu), xl(2:end)];
Z = [flip(zu), zl(2:end)];

N = length(X) - 1;

%% 7. PUNTOS DE CONTROL Y NORMALES
xc = zeros(N,1);
zc_p = zeros(N,1);
nx = zeros(N,1);
nz = zeros(N,1);

% Determinar orientación del contorno
Area = 0.5 * sum( X(1:end-1).*Z(2:end) - X(2:end).*Z(1:end-1) );
isCCW = Area > 0;   % antihorario

for i = 1:N
    dx = X(i+1) - X(i);
    dz = Z(i+1) - Z(i);

    % Punto de control
    xc(i) = 0.5*(X(i) + X(i+1));
    zc_p(i) = 0.5*(Z(i) + Z(i+1));

    L = hypot(dx,dz);

    if isCCW
        nx(i) =  dz / L;
        nz(i) = -dx / L;
    else
        nx(i) = -dz / L;
        nz(i) =  dx / L;
    end
end

%% 8. REPRESENTACIÓN GRÁFICA (EXIGIDA EN A.1)
figure('Color','w','Name','A.1 – Geometría y panelado NACA 6409');
hold on;

plot(X, Z, 'b-o', 'LineWidth',1.2, 'MarkerSize',3);
plot(x_dist, zc, 'r--', 'LineWidth',1.2);
plot(xc, zc_p, 'kx', 'MarkerSize',6);

quiver(xc, zc_p, nx, nz, 0.3, 'Color',[0 0.6 0], 'LineWidth',1.2);

axis equal;
grid on;
xlabel('x/c');
ylabel('z/c');
title(sprintf('Perfil NACA 6409 – Panelado (N = %d)', N));
legend('Contorno','Línea media','Puntos de control','Normales','Location','best');