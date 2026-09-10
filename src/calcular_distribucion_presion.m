%% C.1 – DISTRIBUCION DE PRESION
%% Perfil NACA 6409 – Representación del coeficiente de presión

clc; clear; close all;

%% 1. PARÁMETROS DEL PERFIL (NACA 6409)
f_max = 0.06;   % curvatura máxima (M = 6 → 6%)
x_fmax = 0.40;   % posición curvatura (P = 4 → 40%)
t_max = 0.09;   % espesor máximo (XX = 09 → 9%)
c = 1.0;   % cuerda [m]
U_inf = 1.0;   % velocidad libre [m/s]
rho = 1.225;   % densidad [kg/m³]
n_p = 100;   % paneles por superficie (extradós / intradós)
                   % → N_total = 2*n_p paneles

%% 2. ÁNGULOS DE ATAQUE
alpha_vec_deg = -5:1:10;   % para gráficas de Cl, Cm y Cp

%% 3. A.1 — GENERACIÓN DE LA GEOMETRÍA
% Distribución coseno (concentra nodos en BA y BS)
beta  = linspace(0, pi, n_p + 1);
x_cos = 0.5 * (1 - cos(beta));       % n_p+1 puntos en [0, 1]

% Línea de curvatura media y su derivada
zc_cos    = zeros(size(x_cos));
dzc_cos   = zeros(size(x_cos));
for k = 1:length(x_cos)
    x = x_cos(k);
    if x < x_fmax
        zc_cos(k)  = (f_max/x_fmax^2)     * (2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/x_fmax^2)   * (x_fmax - x);
    else
        zc_cos(k)  = (f_max/(1-x_fmax)^2) * ((1-2*x_fmax) + 2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/(1-x_fmax)^2) * (x_fmax - x);
    end
end

% Distribución de espesor
ze_cos = 5*t_max*(0.2969*sqrt(x_cos) - 0.1260*x_cos - 0.3516*x_cos.^2 ...
                + 0.2843*x_cos.^3 - 0.1015*x_cos.^4);

theta_geo = atan(dzc_cos);

% Extradós e intradós
xu = x_cos - ze_cos .* sin(theta_geo);
zu = zc_cos + ze_cos .* cos(theta_geo);
xl = x_cos + ze_cos .* sin(theta_geo);
zl = zc_cos - ze_cos .* cos(theta_geo);

% Nodos: BS → extradós (↑x=1→x=0) → BA → intradós (↑x=0→x=1) → BS
Xn = [flip(xu),  xl(2:end)];
Zn = [flip(zu),  zl(2:end)];
N  = length(Xn) - 1;    % número total de paneles = 2*n_p

%% 4. Parámetros geométricos de los paneles
xcp  = zeros(N,1);  zcp  = zeros(N,1);
Lp   = zeros(N,1);  phi  = zeros(N,1);
nx_p = zeros(N,1);  nz_p = zeros(N,1);
tx_p = zeros(N,1);  tz_p = zeros(N,1);
for i = 1:N
    dx      = Xn(i+1) - Xn(i);
    dz      = Zn(i+1) - Zn(i);
    Lp(i)   = sqrt(dx^2 + dz^2);
    phi(i)  = atan2(dz, dx);
    xcp(i)  = 0.5*(Xn(i) + Xn(i+1));
    zcp(i)  = 0.5*(Zn(i) + Zn(i+1));
    tx_p(i) =  cos(phi(i));
    tz_p(i) =  sin(phi(i));
    nx_p(i) = -sin(phi(i));
    nz_p(i) =  cos(phi(i));
end

%% 5. A.2 — MATRICES DE INFLUENCIA (método coordenadas locales)

% Notación: An(i,j) influencia NORMAL de FUENTE j en control i
%           Bn(i,j) influencia NORMAL de VÓRTICE j en control i
%           At(i,j) influencia TANGENCIAL de FUENTE j en control i
%           Bt(i,j) influencia TANGENCIAL de VÓRTICE j en control i

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
            dxi  = xcp(i) - Xn(j);
            dzi  = zcp(i) - Zn(j);

            % Sistema local del panel j
            ct = cos(phi(j));   st = sin(phi(j));
            Xi =  ct*dxi + st*dzi;   % coordenada paralela al panel j
            Zi = -st*dxi + ct*dzi;   % coordenada perpendicular

            % Distancias al nodo inicial y final del panel j
            r1sq = Xi^2 + Zi^2;
            r2sq = (Xi - Lp(j))^2 + Zi^2;
            th1 = atan2(Zi,  Xi);
            th2 = atan2(Zi,  Xi - Lp(j));
            dth = th2 - th1;
            lnr = 0.5 * log(r2sq / (r1sq + 1e-30));

            % Componentes en sistema local del panel j

            % Fuente (σ por unidad de longitud)
            u_s_loc =  lnr  / (2*pi);   % paralelo a panel j
            w_s_loc = -dth  / (2*pi);   % perp. a panel j

            % Vórtice (γ constante)
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

%% 6. BUCLE SOBRE ÁNGULOS DE ATAQUE

Cl_pan = zeros(size(alpha_vec_deg));
Cm_pan = zeros(size(alpha_vec_deg));
Cl_tpl = zeros(size(alpha_vec_deg));
Cm_tpl = zeros(size(alpha_vec_deg));
Cp_store = zeros(N, length(alpha_vec_deg));

% TPL: coeficientes de Fourier de dzc/dx
Nt = 1000;
theta_t = linspace(0, pi, Nt+1);
x_t     = 0.5*(1 - cos(theta_t));
dzc_t   = zeros(size(x_t));
for k = 1:length(x_t)
    x = x_t(k);
    if x < x_fmax
        dzc_t(k) = (2*f_max/x_fmax^2)     * (x_fmax - x);
    else
        dzc_t(k) = (2*f_max/(1-x_fmax)^2) * (x_fmax - x);
    end
end
A0_raw = (1/pi) * trapz(theta_t, dzc_t);
A1_tpl = (2/pi) * trapz(theta_t, dzc_t .* cos(theta_t));
A2_tpl = (2/pi) * trapz(theta_t, dzc_t .* cos(2*theta_t));
alpha_L0 = -(A1_tpl/2 - A0_raw);
Cm_tpl_const = (pi/4) * (A2_tpl - A1_tpl);

for k = 1:length(alpha_vec_deg)
    alpha_d = alpha_vec_deg(k);
    alpha_r = deg2rad(alpha_d);

    % Componentes de corriente libre
    Ux = U_inf * cos(alpha_r);
    Uz = U_inf * sin(alpha_r);

    % Construir sistema (N+1) × (N+1)
    M  = zeros(N+1, N+1);
    rh = zeros(N+1, 1);

    % Filas 1..N: condición de no penetración (V·n = 0)
    for i = 1:N
        M(i, 1:N)   = An(i,:);
        M(i, N+1)   = Bn_sum(i);
        rh(i)       = -(Ux * nx_p(i) + Uz * nz_p(i));
    end

    % Fila N+1: condición de Kutta (Vt,1 + Vt,N = 0)
    M(N+1, 1:N)   = At(1,:) + At(N,:);
    M(N+1, N+1)   = Bt_sum(1) + Bt_sum(N);
    rh(N+1)       = -(Ux*(tx_p(1)+tx_p(N)) + Uz*(tz_p(1)+tz_p(N)));

    % Resolver
    sol   = M \ rh;
    sigma = sol(1:N);
    gamma = sol(N+1);

    % Velocidad tangencial y Cp
    Vt = zeros(N,1);
    for i = 1:N
        Vt(i) = Ux*tx_p(i) + Uz*tz_p(i) ...
              + At(i,:)*sigma + gamma*Bt_sum(i);
    end
    Cp = 1 - (Vt/U_inf).^2;
    Cp_store(:,k) = Cp;

    % Cl: teorema de Kutta-Yukowski
    Gamma_tot  = gamma * sum(Lp);
    Cl_pan(k)  = 2 * Gamma_tot / c;

    % Cm_c/4: integración de presiones
    Cm_num = 0;
    for i = 1:N
        Cm_num = Cm_num - Cp(i) * nz_p(i) * (xcp(i) - 0.25*c) * Lp(i);
    end
    Cm_pan(k) = Cm_num / c^2;

    % TPL
    Cl_tpl(k) = 2*pi*(alpha_r - alpha_L0);
    Cm_tpl(k) = Cm_tpl_const;
end

%% 7. GRAFICA — Distribución de Cp vs x/c

figure('Color','white','Position',[50 50 950 520]);
hold on;
cmap = jet(length(alpha_vec_deg));
for k = 1:length(alpha_vec_deg)

    % Extradós: paneles 1..n_p (vienen del flip, van de x=1 a x=0 → ordenar)
    xe = xcp(1:n_p)/c;           Cpe = Cp_store(1:n_p, k);
    xi = xcp(n_p+1:end)/c;       Cpi = Cp_store(n_p+1:end, k);
    [xe_s, ie] = sort(xe);  Cpe_s = Cpe(ie);
    [xi_s, ii] = sort(xi);  Cpi_s = Cpi(ii);
    plot(xe_s, Cpe_s, '-',  'Color', cmap(k,:), 'LineWidth', 1.0);
    plot(xi_s, Cpi_s, '-',  'Color', cmap(k,:), 'LineWidth', 1.0);
end
colormap(jet);
clim([alpha_vec_deg(1), alpha_vec_deg(end)]);

cb = colorbar;
cb.Label.String = 'Ángulo de ataque, \alpha (°)';
cb.Label.FontSize = 11;

cb.Ticks = [-5 0 5 10];
cb.TickLabels = {'-5°','0°','5°','10°'};
set(gca,'YDir','reverse');
xlabel('x/c','FontSize',12); ylabel('C_p','FontSize',12);
title('Distribución de C_p — Método de Paneles, NACA 6409','FontSize',13);
ylim([-12, 2]); xlim([0, 1]);
grid on; box on;