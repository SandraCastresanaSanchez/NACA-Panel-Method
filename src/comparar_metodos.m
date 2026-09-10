%% C.2 - CURVAS POLARES Y ESTUDIO DE CONVERGENCIA
%% Perfil NACA 6409 – Representación de curvas Cl - α, Cm(c/4) - α y Tabla Resumen

clc; clear; close all;

%% 1. PARÁMETROS DEL PERFIL (NACA 6409)
f_max = 0.06;    % curvatura máxima 
x_fmax = 0.40;    % posición curvatura
t_max = 0.09;    % espesor máximo 
c = 1.0;     % cuerda [m]
U_inf = 1.0;     % velocidad libre [m/s]
n_p = 100;     % paneles por superficie (extradós/intradós)

%% 2. ÁNGULOS DE ATAQUE
alpha_vec_deg = -5:1:10;
n_alpha = numel(alpha_vec_deg);

%%  3. A.1 — GEOMETRÍA
beta  = linspace(0, pi, n_p + 1);
x_cos = 0.5 * (1 - cos(beta));

zc_cos  = zeros(size(x_cos));
dzc_cos = zeros(size(x_cos));
for k = 1:length(x_cos)
    x = x_cos(k);
    if x < x_fmax
        zc_cos(k)  = (f_max/x_fmax^2)       * (2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/x_fmax^2)     * (x_fmax - x);
    else
        zc_cos(k)  = (f_max/(1-x_fmax)^2)   * ((1-2*x_fmax) + 2*x_fmax*x - x^2);
        dzc_cos(k) = (2*f_max/(1-x_fmax)^2) * (x_fmax - x);
    end
end

ze_cos = 5*t_max*(0.2969*sqrt(x_cos) - 0.1260*x_cos - 0.3516*x_cos.^2 ...
            + 0.2843*x_cos.^3   - 0.1015*x_cos.^4);

theta_geo = atan(dzc_cos);

xu = x_cos - ze_cos .* sin(theta_geo);
zu = zc_cos + ze_cos .* cos(theta_geo);
xl = x_cos + ze_cos .* sin(theta_geo);
zl = zc_cos - ze_cos .* cos(theta_geo);

% Nodos: BS → extradós (x=1→0) → BA → intradós (0→1) → BS
Xn = [flip(xu), xl(2:end)];
Zn = [flip(zu), zl(2:end)];
N  = length(Xn) - 1;         % total paneles = 2*n_p

%% 4. Parámetros geométricos de panel
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

%% 5. A.2 — MATRICES DE INFLUENCIA
An = zeros(N, N);
Bn = zeros(N, N);
At = zeros(N, N);
Bt = zeros(N, N);

for i = 1:N
    for j = 1:N
        if i == j
            An(i,j) =  0.5;
            At(i,j) =  0.0;
            Bn(i,j) =  0.0;
            Bt(i,j) = -0.5;   % Hess–Smith estándar
        else
            dxi = xcp(i) - Xn(j);
            dzi = zcp(i) - Zn(j);

            ct  = cos(phi(j));  st = sin(phi(j));
            Xi  =  ct*dxi + st*dzi;
            Zi  = -st*dxi + ct*dzi;

            r1sq = Xi^2 + Zi^2;
            r2sq = (Xi - Lp(j))^2 + Zi^2;

            th1 = atan2(Zi, Xi);
            th2 = atan2(Zi, Xi - Lp(j));
            dth = th2 - th1;

            lnr  = 0.5 * log(r2sq / (r1sq + 1e-30));

            % Fuente (σ)
            u_s_loc =  lnr / (2*pi);
            w_s_loc = -dth / (2*pi);

            % Vórtice (γ)
            u_v_loc =  dth / (2*pi);
            w_v_loc =  lnr / (2*pi);

            % Rotar a global
            u_s =  ct*u_s_loc - st*w_s_loc;
            w_s =  st*u_s_loc + ct*w_s_loc;
            u_v =  ct*u_v_loc - st*w_v_loc;
            w_v =  st*u_v_loc + ct*w_v_loc;

            % Proyectar sobre (n,t) del panel i
            An(i,j) = u_s*nx_p(i) + w_s*nz_p(i);
            At(i,j) = u_s*tx_p(i) + w_s*tz_p(i);
            Bn(i,j) = u_v*nx_p(i) + w_v*nz_p(i);
            Bt(i,j) = u_v*tx_p(i) + w_v*tz_p(i);
        end
    end
end

Bn_sum = sum(Bn, 2);
Bt_sum = sum(Bt, 2);

%% 6. TPL — COEFICIENTES DE FOURIER (A0, A1, A2)
Nt      = 1000;
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

A0_raw       = (1/pi) * trapz(theta_t, dzc_t);
A1_tpl       = (2/pi) * trapz(theta_t, dzc_t .* cos(theta_t));
A2_tpl       = (2/pi) * trapz(theta_t, dzc_t .* cos(2*theta_t));

alpha_L0     = -(A1_tpl/2 - A0_raw);                 % rad
Cm_tpl_const = (pi/4) * (A2_tpl - A1_tpl);           % constante (indep. de α)

%% 7. BUCLE ÁNGULOS DE ATAQUE — PANELES + TPL
Cl_pan = zeros(1,n_alpha);
Cm_pan = zeros(1,n_alpha);
Cl_tpl = zeros(1,n_alpha);
Cm_tpl = zeros(1,n_alpha);

for k = 1:n_alpha
    alpha_r = deg2rad(alpha_vec_deg(k));

    Ux = U_inf * cos(alpha_r);
    Uz = U_inf * sin(alpha_r);

    M  = zeros(N+1, N+1);
    rh = zeros(N+1, 1);

    for i = 1:N
        M(i, 1:N) = An(i,:);
        M(i, N+1) = Bn_sum(i);
        rh(i)     = -(Ux*nx_p(i) + Uz*nz_p(i));
    end

    M(N+1, 1:N) = At(1,:) + At(N,:);
    M(N+1, N+1) = Bt_sum(1) + Bt_sum(N);
    rh(N+1)     = -(Ux*(tx_p(1)+tx_p(N)) + Uz*(tz_p(1)+tz_p(N)));

    sol   = M \ rh;
    sigma = sol(1:N);
    gamma = sol(N+1);

    Vt = zeros(N,1);
    for i = 1:N
        Vt(i) = Ux*tx_p(i) + Uz*tz_p(i) + At(i,:)*sigma + gamma*Bt_sum(i);
    end

    Cp = 1 - (Vt/U_inf).^2;

    Cl_pan(k) = 2 * gamma * sum(Lp) / c;

    Cm_num = 0;
    for i = 1:N
        Cm_num = Cm_num - Cp(i)*nz_p(i)*(xcp(i) - 0.25*c)*Lp(i);
    end
    Cm_pan(k) = Cm_num / c^2;

    Cl_tpl(k) = 2*pi*(alpha_r - alpha_L0);
    Cm_tpl(k) = Cm_tpl_const;
end

%% 8. DATOS XFLR5
alpha_xflr5 = [-5,  -4,   -3,   -2,   -1,    0, ...
                1,   2,    3,    4,    5,    6, ...
                7,   8,    9,   10];

Cl_xflr5    = [-0.19, -0.05,  0.09,  0.35,  0.50,  0.63, ...
                0.76,  0.88,  0.98,  1.09,  1.19,  1.28, ...
                1.38,  1.48,  1.50,  1.40];

Cm_xflr5    = [-0.069, -0.092, -0.113, -0.146, -0.152, -0.154, ...
               -0.153, -0.153, -0.152, -0.150, -0.149, -0.147, ...
               -0.144, -0.138, -0.120, -0.096];

%% 9. ERRORES RELATIVOS [%] (TPL y XFLR5 respecto a PANELES)
eCl_tpl = NaN(1,n_alpha);  eCl_xfl = NaN(1,n_alpha);
eCm_tpl = NaN(1,n_alpha);  eCm_xfl = NaN(1,n_alpha);

for k = 1:n_alpha
    if abs(Cl_pan(k)) > 1e-6
        eCl_tpl(k) = abs(Cl_tpl(k)   - Cl_pan(k)) / abs(Cl_pan(k)) * 100;
        eCl_xfl(k) = abs(Cl_xflr5(k) - Cl_pan(k)) / abs(Cl_pan(k)) * 100;
    end
    if abs(Cm_pan(k)) > 1e-6
        eCm_tpl(k) = abs(Cm_tpl(k)   - Cm_pan(k)) / abs(Cm_pan(k)) * 100;
        eCm_xfl(k) = abs(Cm_xflr5(k) - Cm_pan(k)) / abs(Cm_pan(k)) * 100;
    end
end

%% 10. GRAFICA 1 - COMPARATIVA C_l Y C_m vs α
figure('Color','white','Position',[50 50 1100 460]);

subplot(1,2,1)
plot(alpha_vec_deg,Cl_pan,'b-o','LineWidth',2,'MarkerFaceColor','b'); hold on
plot(alpha_vec_deg,Cl_tpl,'k--','LineWidth',2)
grid on; box on
xlim([-5 10]); ylim([0 2.5])
xlabel('α (deg)'); ylabel('C_l')
title('Comparativa C_l vs \alpha — NACA 6409')

subplot(1,2,2)
plot(alpha_vec_deg,Cm_pan,'g-s','LineWidth',2,'MarkerFaceColor','g'); hold on
plot(alpha_vec_deg,Cm_tpl,'k--','LineWidth',2)
grid on; box on
xlim([-5 10]); ylim([-0.22 -0.12])
xlabel('α (deg)'); ylabel('C_{m c/4}')
title('Comparativa C_m vs \alpha — NACA 6409')

%% 11. GRAFICA 2 — Cl vs α
figure('Color','white','Position',[100 100 820 560]);
plot(alpha_vec_deg, Cl_pan, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Método Paneles (Hess-Smith)');
hold on;
plot(alpha_vec_deg, Cl_tpl, 'k--', 'LineWidth', 2, ...
     'DisplayName', 'TPL');
plot(alpha_xflr5, Cl_xflr5, 'r-^', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'XFLR5 (viscoso)');
yline(0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 0.8, ...
     'HandleVisibility', 'off');
xlabel('\alpha (deg)', 'FontSize', 13);
ylabel('C_l',          'FontSize', 13);
title('Comparativa C_l vs \alpha — NACA 6409', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
xlim([-5, 10]); ylim([-0.4, 2.2]);
grid on; box on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

%% 12. GRAFICA 3 — Cm(c/4) vs α
figure('Color','white','Position',[120 120 820 560]);
plot(alpha_vec_deg, Cm_pan, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Método Paneles (Hess-Smith)');
hold on;
plot(alpha_vec_deg, Cm_tpl, 'k--', 'LineWidth', 2, ...
     'DisplayName', 'TPL');
plot(alpha_xflr5, Cm_xflr5, 'r-^', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'XFLR5 (viscoso)');
yline(0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 0.8, ...
     'HandleVisibility', 'off');
xlabel('\alpha (deg)', 'FontSize', 13);
ylabel('C_{m_{c/4}}',  'FontSize', 13);
title('Comparativa C_m vs \alpha — NACA 6409', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
xlim([-5, 10]); ylim([-0.20, 0.02]);
grid on; box on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

%% 13. TABLA VALORES
fmtNA = @(v,fmt) ternary(isnan(v), 'N/A', sprintf(fmt, v));

col_headers = {'\alpha (°)', ...
               'Cl Paneles','Cl TPL','Cl XFLR5', ...
               'Err Cl TPL (%)','Err Cl XFLR5 (%)', ...
               'Cm Paneles','Cm TPL','Cm XFLR5', ...
               'Err Cm TPL (%)','Err Cm XFLR5 (%)'};

data_cell = cell(n_alpha, 11);

for k = 1:n_alpha
    data_cell{k,1}  = sprintf('%d', alpha_vec_deg(k));  
    data_cell{k,2}  = sprintf('%.4f', Cl_pan(k));       
    data_cell{k,3}  = sprintf('%.4f', Cl_tpl(k));  
    data_cell{k,4}  = sprintf('%.2f', Cl_xflr5(k));     
    data_cell{k,5}  = fmtNA(eCl_tpl(k), '%.2f');        
    data_cell{k,6}  = fmtNA(eCl_xfl(k), '%.2f');
    data_cell{k,7}  = sprintf('%.4f', Cm_pan(k));       
    data_cell{k,8}  = sprintf('%.4f', Cm_tpl(k));       
    data_cell{k,9}  = sprintf('%.3f', Cm_xflr5(k));    
    data_cell{k,10} = fmtNA(eCm_tpl(k), '%.2f');        
    data_cell{k,11} = fmtNA(eCm_xfl(k), '%.2f');              
end

figure('Color','white','Position',[30 30 1430 570], 'Name','Tabla Resumen NACA 6409');
axis off;

uitable('Data', data_cell, 'ColumnName', col_headers, ...
        'RowName', [], 'Units', 'normalized', 'Position', ...
        [0.005 0.04 0.99 0.87], 'FontSize', 10, 'ColumnWidth', ...
        {38,82,68,68,108,118,82,68,72,108,118});

annotation('textbox',[0 0.92 1 0.08], ...
    'String', ['Tabla Resumen — NACA 6409: C_l y C_{m c/4} (Paneles | TPL | XFLR5) ' ...
               '— Errores relativos respecto al Método de Paneles'], ...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','EdgeColor','none');

%% helper ternary
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

%% 14. PRINTS
fprintf('alpha_L0 (TPL) = %.3f deg\n', rad2deg(alpha_L0));
fprintf('Cm_tpl (constante) = %.4f\n', Cm_tpl_const);
fprintf('Cl_max XFLR5 = %.2f en alpha = %d deg\n', max(Cl_xflr5), alpha_xflr5(Cl_xflr5 == max(Cl_xflr5)));
fprintf('Cm_xflr5 min = %.3f en alpha = %d deg\n', min(Cm_xflr5), alpha_xflr5(Cm_xflr5 == min(Cm_xflr5)));