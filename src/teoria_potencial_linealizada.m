%% B – TEORIA POTENCIAL LINEALIZADA DE PERFILES EN REGIMEN INCOMPRESIBLE
%% Perfil NACA 6409 – Calcular A0, A1, A2 y tabla Cl, Cm(c/4) vs α

clear; clc; close all;

%% 1. PARÁMETROS DEL PERFIL (NACA 6409)
f_max = 0.06;   % curvatura máxima
x_fmax = 0.40;   % posición de curvatura máxima

%% 2. RESOLUCION DE LA INTEGRAL EN θ
nTheta = 5000;

%% 3. FUNCIÓN: coeficientes A0, A1, A2 por integración (TPL)
function [A0_camber, A1, A2] = TPL_Fourier_NACA4_camber(f_max, x_fmax, nTheta)

% Integra en θ ∈ [0,π] usando x = (1 - cosθ)/2
% dzc/dx definido a trozos para NACA 4 cifras (línea media)

    if nargin < 3, nTheta = 4000; end
    theta = linspace(0, pi, nTheta);
    x = 0.5*(1 - cos(theta));   % mapeo TPL

    % dzc/dx a trozos
    dzdx = zeros(size(x));
    i1 = (x < x_fmax);
    i2 = ~i1;

    dzdx(i1) = 2*f_max/(x_fmax^2) * (x_fmax - x(i1));
    dzdx(i2) = 2*f_max/((1-x_fmax)^2) * (x_fmax - x(i2));

    % Coeficientes de Fourier
    A0_camber = -(1/pi) * trapz(theta, dzdx);
    A1 = (2/pi) * trapz(theta, dzdx .* cos(theta));
    A2 = (2/pi) * trapz(theta, dzdx .* cos(2*theta));
end

%% 4. COEFICIENTES DE FOURIER de TPL
[A0_camber, A1, A2] = TPL_Fourier_NACA4_camber(f_max, x_fmax, nTheta);

fprintf("A0_camber = %.6f\nA1 = %.6f\nA2 = %.6f\n\n", A0_camber, A1, A2);

%% 5. BARRIDO DE ANGULOS DE ATAQUE
alpha_deg = (-5:10).';
alpha_rad = deg2rad(alpha_deg);

A0 = alpha_rad + A0_camber;

Cl_TPL = pi*(2*A0 + A1);
Cm_c4_TPL = (pi/4)*(A2 - A1) * ones(size(alpha_deg));

%% 6. TABLA GRÁFICA

n = length(alpha_deg);

col_headers = {'α (°)', 'A₀', 'C_l', 'C_{m,c/4}'};
data_cell = cell(n, 4);

for k = 1:n
    data_cell{k,1} = sprintf('%d', alpha_deg(k));
    data_cell{k,2} = sprintf('%.5f', A0(k));
    data_cell{k,3} = sprintf('%.5f', Cl_TPL(k));
    data_cell{k,4} = sprintf('%.5f', Cm_c4_TPL(k));
end

figure('Color','white','Position',[50 50 720 420], ...
       'Name','Tabla TPL — NACA 6409');
axis off;

uitable('Data', data_cell, 'ColumnName', col_headers, ...
        'RowName', [], 'Units', 'normalized', 'Position', ...
        [0.04 0.08 0.92 0.78], 'FontSize', 11, ...
        'ColumnWidth', {80, 130, 130, 150});

annotation('textbox',[0 0.96 1 0.04], ...
    'String','TPL — NACA 6409', ...
    'HorizontalAlignment','center', ...
    'FontSize',14,'FontWeight','bold', ...
    'EdgeColor','none');

coef_text = sprintf(['Coeficientes TPL del perfil NACA 6409:\n' ...
                     'A_{0,camber} = %.5f\n' ...
                     'A_1 = %.5f\n' ...
                     'A_2 = %.5f'], ...
                     A0_camber, A1, A2);

annotation('textbox',[0.32 0.88 0.36 0.08], ...
    'String', coef_text, ...
    'FitBoxToText','on', ...
    'BackgroundColor','white', ...
    'FontSize',10);

fprintf('Tabla gráfica TPL generada correctamente.\n');