%% A.2 – RESOLUCION DEL SISTEMA DE ECUACIONES 
%% Perfil NACA 6409 – Fuentes σ y vorticidad γ

clear; clc; close all;

%% 1. PARÁMETROS DEL PERFIL (NACA 6409)
f_max = 0.06;   % curvatura máxima
x_fmax = 0.40;   % posición curvatura
t_max = 0.09;   % espesor máximo
c = 1.0;   % cuerda [m]
U_inf = 1.0;   % velocidad libre [m/s]
rho = 1.225;   % densidad [kg/m³]
n_p = 100;   % paneles por superficie (extradós / intradós)
             % → N_total = 2*n_p paneles

%% 2. ÁNGULO DE ATAQUE
alpha_deg = 5;

%% 3. A.1 — GENERACIÓN DE LA GEOMETRÍA

% Distribución coseno (concentra nodos en BA y BS)
beta  = linspace(0, pi, n_p + 1);
x_cos = 0.5 * (1 - cos(beta));

% Línea de curvatura media y su derivada
zc_cos  = zeros(size(x_cos));
dzc_cos = zeros(size(x_cos));
for k = 1:length(x_cos)
    x = x_cos(k);
    if x < x_fmax
        zc_cos(k) = (f_max/x_fmax^2) * (2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/x_fmax^2) * (x_fmax - x);
    else
        zc_cos(k) = (f_max/(1-x_fmax)^2) * ((1-2*x_fmax) + 2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/(1-x_fmax)^2) * (x_fmax - x);
    end
end

% Distribución de espesor
ze_cos = 5*t_max*(0.2969*sqrt(x_cos) - 0.1260*x_cos - 0.3516*x_cos.^2 ...
                + 0.2843*x_cos.^3 - 0.1036*x_cos.^4);

theta_geo = atan(dzc_cos);

% Extradós e intradós
xu = x_cos - ze_cos .* sin(theta_geo);
zu = zc_cos + ze_cos .* cos(theta_geo);
xl = x_cos + ze_cos .* sin(theta_geo);
zl = zc_cos - ze_cos .* cos(theta_geo);

% Nodos: BS → extradós (↑x=1→x=0) → BA → intradós (↑x=0→x=1) → BS
Xn = [flip(xu),  xl(2:end)];
Zn = [flip(zu),  zl(2:end)];
N  = length(Xn) - 1;   % número total de paneles = 2*n_p

%% 4. Parámetros geométricos de los paneles
xcp  = zeros(N,1);  zcp = zeros(N,1);
Lp   = zeros(N,1);  phi = zeros(N,1);
nx_p = zeros(N,1);  nz_p = zeros(N,1);
tx_p = zeros(N,1);  tz_p = zeros(N,1);

for i = 1:N
    dx = Xn(i+1) - Xn(i);
    dz = Zn(i+1) - Zn(i);
    Lp(i) = sqrt(dx^2 + dz^2);
    phi(i) = atan2(dz, dx);
    xcp(i) = 0.5*(Xn(i) + Xn(i+1));
    zcp(i) = 0.5*(Zn(i) + Zn(i+1));
    tx_p(i) = cos(phi(i));
    tz_p(i) = sin(phi(i));
    nx_p(i) = -sin(phi(i));
    nz_p(i) = cos(phi(i));
end

%% 5. A.2 — MATRICES DE INFLUENCIA (método coordenadas locales)
An = zeros(N, N);
Bn = zeros(N, N);
At = zeros(N, N);
Bt = zeros(N, N);

for i = 1:N
    for j = 1:N
        if i == j
            % Singularidad: límite analítico conocido
            An(i,j) =  0.5;
            At(i,j) =  0.0;
            Bn(i,j) =  0.0;
            Bt(i,j) = -0.5;   % signo estándar Hess-Smith
        else
            % Vector del nodo j al punto de control i en sistema global
            dxi = xcp(i) - Xn(j);
            dzi = zcp(i) - Zn(j);

            % Sistema local del panel j
            ct = cos(phi(j));   st = sin(phi(j));
            Xi = ct*dxi + st*dzi;   % coordenada paralela al panel j
            Zi = -st*dxi + ct*dzi;   % coordenada perpendicular

            % Distancias al nodo inicial y final del panel j
            r1sq = Xi^2 + Zi^2;
            r2sq = (Xi - Lp(j))^2 + Zi^2;

            th1 = atan2(Zi,  Xi);
            th2 = atan2(Zi,  Xi - Lp(j));
            dth = th2 - th1;

            lnr = 0.5 * log(r2sq / (r1sq + 1e-30));

            % Componentes en sistema local del panel j
            % — Fuente (σ por unidad de longitud) —
            u_s_loc =  lnr  / (2*pi);   % paralelo a panel j
            w_s_loc = -dth  / (2*pi);   % perp. a panel j

            % — Vórtice (γ constante) —
            u_v_loc =  dth  / (2*pi);
            w_v_loc =  lnr  / (2*pi);

            % Rotar al sistema global
            u_s =  ct*u_s_loc - st*w_s_loc;
            w_s =  st*u_s_loc + ct*w_s_loc;
            u_v =  ct*u_v_loc - st*w_v_loc;
            w_v =  st*u_v_loc + ct*w_v_loc;

            % Proyectar sobre normal y tangente del panel i
            An(i,j) =  u_s * nx_p(i) + w_s * nz_p(i);
            At(i,j) =  u_s * tx_p(i) + w_s * tz_p(i);
            Bn(i,j) =  u_v * nx_p(i) + w_v * nz_p(i);
            Bt(i,j) =  u_v * tx_p(i) + w_v * tz_p(i);
        end
    end
end

% Vectores suma de columnas (coeficiente de γ único)
Bn_sum = sum(Bn, 2);   % (N×1)
Bt_sum = sum(Bt, 2);   % (N×1)

%% 6. RESOLVER PARA EL ÁNGULO ESPECIFICADO
alpha_r = deg2rad(alpha_deg);

% Componentes de corriente libre
Ux = U_inf * cos(alpha_r);
Uz = U_inf * sin(alpha_r);

% Construir sistema (N+1) × (N+1)
M  = zeros(N+1, N+1);
rh = zeros(N+1, 1);

% Filas 1..N: condición de no penetración (V·n = 0)
for i = 1:N
    M(i, 1:N) = An(i,:);
    M(i, N+1) = Bn_sum(i);
    rh(i) = -(Ux * nx_p(i) + Uz * nz_p(i));
end

% Fila N+1: condición de Kutta (Vt,1 + Vt,N = 0)
M(N+1, 1:N) = At(1,:) + At(N,:);
M(N+1, N+1) = Bt_sum(1) + Bt_sum(N);
rh(N+1) = -(Ux*(tx_p(1)+tx_p(N)) + Uz*(tz_p(1)+tz_p(N)));

% Resolver
sol = M \ rh;
sigma = sol(1:N);
gamma = sol(N+1);

%% 7. Gráfica — Distribución de fuentes σ y vorticidad γ
% Extradós (paneles 1..n_p) e intradós (n_p+1..N): ordenar por x/c creciente
% Se muestra -sigma por la misma convención de signo que gamma
xe  = xcp(1:n_p)/c;        sig_e = -sigma(1:n_p);
xi  = xcp(n_p+1:end)/c;    sig_i = -sigma(n_p+1:end);
[xe_s, ie] = sort(xe);   sig_e_s = sig_e(ie);
[xi_s, ii] = sort(xi);   sig_i_s = sig_i(ii);

figure('Color','white','Position',[50 50 820 560]);
hold on;

% Fuentes extradós — línea continua azul
plot(xe_s, sig_e_s, 'b-',  'LineWidth', 2.0, 'DisplayName', 'Fuentes q (Extradós)');

% Fuentes intradós — línea discontinua azul
plot(xi_s, sig_i_s, 'b--', 'LineWidth', 2.0, 'DisplayName', 'Fuentes q (Intradós)');

% Vorticidad γ = constante — línea roja horizontal
% Se muestra -gamma porque con esta convención gamma<0 para sustentación positiva
gamma_display = -gamma;
yline(gamma_display, 'r-', 'LineWidth', 1.8, 'Label', ...
      sprintf('  Vorticidad γ = %.4f', gamma_display), ...
      'LabelVerticalAlignment','bottom', 'DisplayName', ...
      sprintf('Vorticidad γ = %.4f', gamma_display));

xlabel('x/c','FontSize',12);
ylabel('Intensidad de Singularidad','FontSize',12);
title(sprintf('Distribución de Fuentes y Vorticidad (α=%d°)', alpha_deg),'FontSize',13);
legend('Location','northwest','FontSize',10);
grid on; box on;
xlim([0, 1]);

fprintf('\n--- Resultados a α = %d° ---\n', alpha_deg);
fprintf('γ = %.4f\n', -gamma);
fprintf('\n--- Programa finalizado ---\n');