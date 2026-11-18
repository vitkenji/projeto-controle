%% ========================================================================
%  UTFPR - DAELN | Controle 1 / Laboratório de Controle (2026/2)
%  Alunos: Vitor Kenji Zoppo Yamada e Adryan Castro Feres
%  Script completo (itens a até f) alinhado à notação dos slides da
%  disciplina, com cálculos idênticos e gráficos em alta definição.
%% ========================================================================
clear; clc; close all;

% Garante tema claro (Light Mode) como padrão global antes de criar figuras
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');

%% PARÂMETROS EXCLUSIVOS DA EQUIPE: Vitor Kenji Zoppo Yamada e Adryan Castro Feres
K_planta  = 160;
a         = 8;
b         = 20;

Kv_des    = 10;     % Constante estática de erro de velocidade desejada (s^-1)
MG_des_dB = 30;     % Margem de Ganho mínima desejada (dB)
MF_des    = 80;     % Margem de Fase mínima desejada (graus)

fprintf('=========================================================================\n');
fprintf('PARÂMETROS DA EQUIPE: Vitor Kenji Zoppo Yamada e Adryan Castro Feres\n');
fprintf('K = %g | a = %g | b = %g | Kv_des = %g s^-1 | MG_des = %g dB | MF_des = %g deg\n', ...
    K_planta, a, b, Kv_des, MG_des_dB, MF_des);
fprintf('=========================================================================\n\n');

%% ========================================================================
%  ITEM (a): DIAGRAMA DE LUGAR DAS RAÍZES, BODE E ANÁLISE DE ESTABILIDADE
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (a): DIAGRAMAS DE LGR E BODE E ANÁLISE DE ESTABILIDADE\n');
fprintf('*************************************************************************\n');

s = tf('s');
num_G = K_planta;
den_G = conv([1 0], conv([1 a], [1 b])); % s*(s+8)*(s+20) = s^3 + 28s^2 + 160s
G = tf(num_G, den_G);

polos_G = pole(G);
zeros_G = zero(G);
n_polos = length(polos_G);
m_zeros = length(zeros_G);
Kv_orig = K_planta / (a * b);

% Cálculos analíticos do Lugar das Raízes (LGR)
centroide_LGR = (sum(polos_G) - sum(zeros_G)) / (n_polos - m_zeros);
angulos_assintotas = ((2*(0:(n_polos - m_zeros - 1)) + 1) * 180) / (n_polos - m_zeros);

% Ponto de quebra (breakaway): d/ds [s^3 + 28s^2 + 160s] = 3s^2 + 56s + 160 = 0
deriv_den = polyder(den_G);
raizes_break = roots(deriv_den);
s_break_valido = raizes_break(raizes_break < 0 & raizes_break > -a);
K_cr_break = -polyval(den_G, s_break_valido) / K_planta;

% Valores críticos via Critério de Routh-Hurwitz:
% s^3 + 28s^2 + 160s + 160*K_cr = 0 => K_cr = (28 * 160) / 160 = 28
K_cr = (den_G(2) * den_G(3)) / K_planta;
K_total_crit = K_planta * K_cr;
w_crit_imag = sqrt(den_G(3)); % sqrt(160) = 12.6491 rad/s

% Margens de estabilidade originais no Bode (w_cf = cruz. fase, w_cg = cruz. ganho)
[MG_orig_abs, MF_orig, w_cf_orig, w_cg_orig] = margin(G);
MG_orig_dB = 20 * log10(MG_orig_abs);
[mag_wcf_orig, phase_wcf_orig] = bode(G, w_cf_orig);
[mag_wcg_orig, phase_wcg_orig] = bode(G, w_cg_orig);
mag_wcf_orig_dB = 20 * log10(mag_wcf_orig);

fprintf('1. Função de Transferência G(s):\n');
fprintf('   Numerador   : %g\n', num_G);
fprintf('   Denominador : %.0f*s^3 + %.0f*s^2 + %.0f*s + %.0f\n', den_G(1), den_G(2), den_G(3), den_G(4));
fprintf('   Polos de MA : %.4f, %.4f, %.4f\n', polos_G(1), polos_G(2), polos_G(3));
fprintf('   Kv original : %.4f s^-1 (erro à rampa original ess = %.4f)\n\n', Kv_orig, 1/Kv_orig);

fprintf('2. Parâmetros Intermediários do Lugar das Raízes (LGR):\n');
fprintf('   Centroide das assíntotas (sigma_a) : %.4f\n', centroide_LGR);
fprintf('   Ângulos das assíntotas             : %.1f deg, %.1f deg, %.1f deg\n', angulos_assintotas);
fprintf('   Raízes de d(den)/ds = 0            : %.4f e %.4f\n', raizes_break(1), raizes_break(2));
fprintf('   Ponto de quebra válido no LGR      : s = %.4f (com ganho K = %.4f)\n\n', ...
    s_break_valido, K_cr_break);

fprintf('3. Valores Críticos de Ganho, Fase e Frequência:\n');
fprintf('   Ganho crítico adicional (K_cr)         : %.4f\n', K_cr);
fprintf('   Ganho total crítico (K_total_crit)     : %.4f\n', K_total_crit);
fprintf('   Freq. cruzamento de fase (w_cf)        : %.4f rad/s (polos em +-j%.4f)\n', w_cf_orig, w_crit_imag);
fprintf('   Módulo |G(j*w_cf)|                     : %.6f (%.4f dB) | Fase = %.2f deg\n', ...
    mag_wcf_orig, mag_wcf_orig_dB, phase_wcf_orig);
fprintf('   Margem de Ganho original (MG_orig)     : %.4f V/V (%.4f dB)\n', MG_orig_abs, MG_orig_dB);
fprintf('   Freq. cruzamento de ganho (w_cg)       : %.4f rad/s | |G(j*w_cg)| = %.4f (0 dB)\n', w_cg_orig, mag_wcg_orig);
fprintf('   Fase em w_cg (<G(j*w_cg))              : %.4f deg\n', phase_wcg_orig);
fprintf('   Margem de Fase original (MF_orig)      : %.4f deg\n\n', MF_orig);

fprintf('4. Análise de Estabilidade (Item a):\n');
fprintf('   O sistema em malha fechada unitária sem compensação é ESTÁVEL, pois MG_orig = %.2f dB > 0 dB\n', MG_orig_dB);
fprintf('   (0 < K = 1 < %.2f), MF_orig = %.2f deg > 0 deg e w_cg (%.2f rad/s) < w_cf (%.2f rad/s).\n', ...
    K_cr, MF_orig, w_cg_orig, w_cf_orig);
fprintf('   Contudo, Kv_orig = %.2f s^-1 não atende à especificação de Kv = %g s^-1.\n\n', Kv_orig, Kv_des);

% -------------------------------------------------------------------------
% Gráfico (a.1): Lugar das Raízes (padrão visual idêntico ao modelo)
% -------------------------------------------------------------------------
fig_a1 = figure('Name', 'Item (a) - Lugar das Raízes', 'Position', [100 100 840 660]);
k_lgr_a = unique([0, logspace(-4, 3.5, 4000), K_cr_break, K_cr]);
r_lgr_a = rlocus(G, k_lgr_a);

hold on;
% Eixos cartesianos de referência
xline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');

% Ramos do LGR
plot(real(r_lgr_a(1,:)), imag(r_lgr_a(1,:)), '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2, 'DisplayName', 'Ramos LGR');
for idx_r = 2:size(r_lgr_a, 1)
    plot(real(r_lgr_a(idx_r,:)), imag(r_lgr_a(idx_r,:)), '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2, 'HandleVisibility', 'off');
end

% Polos de Malha Aberta, Centroide, Ponto de Quebra e Cruzamento jw
plot(real(polos_G), imag(polos_G), 'rx', 'MarkerSize', 9, 'LineWidth', 2.2, 'DisplayName', 'Polos MA');
plot(centroide_LGR, 0, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Centroide');
plot(s_break_valido, 0, 'o', 'Color', [0 0.55 0.2], 'MarkerFaceColor', [0 0.55 0.2], 'MarkerSize', 8, 'DisplayName', 'Ponto de Quebra');
plot([0 0], [w_crit_imag -w_crit_imag], 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8, 'DisplayName', 'Cruzamento Eixo j\omega');

% Assíntotas a partir do centroide
r_assint = 35;
for idx_ang = 1:length(angulos_assintotas)
    ang_rad = deg2rad(angulos_assintotas(idx_ang));
    x_ass = [centroide_LGR, centroide_LGR + r_assint * cos(ang_rad)];
    y_ass = [0, r_assint * sin(ang_rad)];
    if idx_ang == 1
        plot(x_ass, y_ass, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, 'DisplayName', 'Assíntotas');
    else
        plot(x_ass, y_ass, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end
end

% Rótulos diretos dos pontos notáveis
text(centroide_LGR, -1.15, sprintf('Centroide (%.2f)', centroide_LGR), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k', 'FontSize', 10);
text(s_break_valido, 1.15, sprintf('Quebra (%.2f, K=%.2f)', s_break_valido, K_cr_break), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0 0.55 0.2], 'FontSize', 10);
text(0.45, w_crit_imag, sprintf('j%.2f (K_{cr}=%.0f)', w_crit_imag, K_cr), ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'Color', [0.55 0 0.65], 'FontSize', 10);
text(0.45, -w_crit_imag, sprintf('-j%.2f (K_{cr}=%.0f)', w_crit_imag, K_cr), ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'Color', [0.55 0 0.65], 'FontSize', 10);

xlim([-25 5]); ylim([-18 18]);
xticks(-25:5:5); yticks(-15:5:15);
title('Lugar das Raízes de G(s) = 160 / [s(s+8)(s+20)]', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Eixo Real (s^{-1})', 'FontSize', 11);
ylabel('Eixo Imaginário (rad/s)', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9.5);
grid on; grid minor; box off;
aplicar_light_mode_e_exportar(fig_a1, 'item_a_LGR.png');

% -------------------------------------------------------------------------
% Gráfico (a.2): Diagrama de Bode (padrão visual idêntico ao modelo)
% -------------------------------------------------------------------------
fig_a2 = figure('Name', 'Item (a) - Diagrama de Bode', 'Position', [140 120 860 680]);
w_bode_a = logspace(-1, 2.4, 2500);
[mag_a_raw, phase_a_raw] = bode(G, w_bode_a);
mag_a_dB = 20 * log10(squeeze(mag_a_raw));
phase_a_deg = unwrap(squeeze(phase_a_raw) * pi/180) * 180/pi;
if phase_a_deg(1) > 0, phase_a_deg = phase_a_deg - 360; end

% Subgráfico Superior: Magnitude (dB)
ax_a2_mag = subplot(2, 1, 1);
semilogx(w_bode_a, mag_a_dB, '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2); hold on;
yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
plot(w_cg_orig, 0, 'o', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 8);
text(w_cg_orig * 1.18, 9, sprintf('\\omega_{cg} = %.2f rad/s (0 dB)', w_cg_orig), ...
    'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 10);

plot([w_cf_orig w_cf_orig], [mag_wcf_orig_dB 0], '-', 'Color', [0.55 0 0.65], 'LineWidth', 2.2);
plot(w_cf_orig, mag_wcf_orig_dB, 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8);
text(w_cf_orig * 1.18, mag_wcf_orig_dB / 2, sprintf('MG = %.2f dB', MG_orig_dB), ...
    'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 10);

xlim([0.1 250]); ylim([-100 40]);
xticks([1e-1 1e0 1e1 1e2]); yticks(-100:20:40);
ylabel('Magnitude (dB)', 'FontSize', 11);
title('Diagrama de Bode de G(s) = 160 / [s(s+8)(s+20)]', 'FontSize', 13, 'FontWeight', 'bold');
grid on; grid minor; box on;

% Subgráfico Inferior: Fase (graus)
ax_a2_pha = subplot(2, 1, 2);
semilogx(w_bode_a, phase_a_deg, '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2); hold on;
yline(-180, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
plot([w_cg_orig w_cg_orig], [-180 phase_wcg_orig], '-', 'Color', [0.8 0 0], 'LineWidth', 2.2);
plot(w_cg_orig, phase_wcg_orig, 'o', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 8);
text(w_cg_orig * 1.18, (-180 + phase_wcg_orig) / 2, sprintf('MF = %.2f^{\\circ}', MF_orig), ...
    'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 10);

plot(w_cf_orig, -180, 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8);
text(w_cf_orig * 1.18, -170, sprintf('\\omega_{cf} = %.2f rad/s (-180^{\\circ})', w_cf_orig), ...
    'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 10);

xlim([0.1 250]); ylim([-278 -80]);
xticks([1e-1 1e0 1e1 1e2]); yticks(-270:45:-90);
xlabel('Frequência (rad/s)', 'FontSize', 11);
ylabel('Fase (graus)', 'FontSize', 11);
grid on; grid minor; box on;
linkaxes([ax_a2_mag, ax_a2_pha], 'x');
aplicar_light_mode_e_exportar(fig_a2, 'item_a_Bode.png');


%% ========================================================================
%  ITEM (b): PROJETO NO DOMÍNIO DA FREQUÊNCIA (AVANÇO-ATRASO - AULA 15)
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (b): CONTROLADOR NO DOMÍNIO DA FREQUÊNCIA D(s) = D_AT(s) * D_AV(s)\n');
fprintf('*************************************************************************\n');

% Passo 1: Avaliar planta apenas com ganho ajustado para Kv = 10 (Kd0 = 10)
Kd_Kv = Kv_des / Kv_orig;
G_Kv = Kd_Kv * G;
[MG_Kv_abs, MF_Kv, wcf_Kv, wcg_Kv] = margin(G_Kv);
fprintf('1. Verificação apenas com ganho proporcional Kd0 = %.2f para Kv = %g:\n', Kd_Kv, Kv_des);
fprintf('   MG resultante = %.4f dB (w_cf = %.4f rad/s) -> Reprova (exigido >= %g dB)\n', ...
    20*log10(MG_Kv_abs), wcf_Kv, MG_des_dB);
fprintf('   MF resultante = %.4f deg (w_cg = %.4f rad/s) -> Reprova (exigido >= %g deg)\n\n', ...
    MF_Kv, wcg_Kv, MF_des);

% Passo 2: Projeto do Compensador Avanço-Atraso em Cascata (Aula 15, Slides 19-23)
% Forma polo-zero: D(s) = Kc_b * [(s + z_av)/(s + p_av)] * [(s + z_at)/(s + p_at)]
% Forma const. tempo: D(s) = Kd * [(tau_av*s + 1)/(alpha*tau_av*s + 1)] * [(tau_at*s + 1)/(alpha_at*tau_at*s + 1)]
z_av = 1.5;
p_av = 2.5;
alpha_av = z_av / p_av;          % alpha < 1 (avanço de fase)
tau_av = 1 / z_av;
phi_max_av = asind((1 - alpha_av) / (1 + alpha_av));
wm_av = sqrt(z_av * p_av);

G_com_av = G * tf([1 z_av], [1 p_av]);
[~, ~, wcf_av, ~] = margin(G_com_av);
[mag_wcf_av, ~] = bode(G_com_av, wcf_av);

% Alvo de MG em alta frequência com folga de +0.8 dB para absorver o pequeno
% atraso de fase (~ -0.5 deg) introduzido pela rede de atraso em w_cf
MG_alvo_projeto_dB = 30.80;
mag_wcf_desejado = 10^(-MG_alvo_projeto_dB / 20);
Kc_b = mag_wcf_desejado / mag_wcf_av;

% Rede de atraso de fase D_AT(s): eleva Kv até 10 sem degradar MG nem MF
alpha_at = Kv_des / (Kc_b * alpha_av * Kv_orig); % alpha_at > 1 (atraso de fase)
z_at = 0.10;
p_at = z_at / alpha_at;
tau_at = 1 / z_at;
Kd_b = Kc_b * alpha_av * alpha_at; % Ganho estático de baixa frequência (Kd = Kv / Kv_orig = 12.5)

D_AV = tf([1 z_av], [1 p_av]);
D_AT = tf([1 z_at], [1 p_at]);
D_b  = Kc_b * D_AV * D_AT;
L_b  = series(D_b, G);
T_b  = feedback(L_b, 1);

% Verificação dinâmica das especificações do item (b)
Kv_b = dcgain(s * L_b);
[MG_b_abs, MF_b, wcf_b, wcg_b] = margin(L_b);
MG_b_dB = 20 * log10(MG_b_abs);
[mag_wcf_b, ~] = bode(L_b, wcf_b);
[~, phase_wcg_b] = bode(L_b, wcg_b);
mag_wcf_b_dB = 20 * log10(mag_wcf_b);
polos_Tb = pole(T_b);
info_b = stepinfo(T_b);

status_Kv_b = verificar_status(abs(Kv_b - Kv_des) < 1e-3);
status_MG_b = verificar_status(MG_b_dB >= MG_des_dB);
status_MF_b = verificar_status(MF_b >= MF_des);

fprintf('2. Parâmetros Intermediários do Compensador Avanço-Atraso D(s) = D_AT(s)*D_AV(s):\n');
fprintf('   Seção de Avanço D_AV(s) : z_av = %.4f | p_av = %.4f | alpha = %.4f (< 1) | tau_av = %.4f s\n', ...
    z_av, p_av, alpha_av, tau_av);
fprintf('   Avanço máximo de fase   : phi_max = %.4f deg em w_max = %.4f rad/s\n', phi_max_av, wm_av);
fprintf('   Freq. w_cf (com D_AV)   : %.4f rad/s | |G*D_AV(jw_cf)| = %.6f\n', wcf_av, mag_wcf_av);
fprintf('   MG alvo (com folga)     : %.2f dB -> |L(jw_cf)|_alvo = %.6f\n', MG_alvo_projeto_dB, mag_wcf_desejado);
fprintf('   Ganho alta freq. (Kc_b) : %.6f | Ganho baixa freq. (Kd) = %.4f\n', Kc_b, Kd_b);
fprintf('   Seção de Atraso D_AT(s) : z_at = %.4f | p_at = %.6f | alpha_at = %.4f (> 1) | tau_at = %.4f s\n\n', ...
    z_at, p_at, alpha_at, tau_at);

[num_Db, den_Db] = tfdata(D_b, 'v');
fprintf('3. Função de Transferência do Controlador D(s):\n');
fprintf('   Forma Const. Tempo : D(s) = %.2f * [(%.4f*s + 1)/(%.4f*s + 1)] * [(%.2f*s + 1)/(%.2f*s + 1)]\n', ...
    Kd_b, tau_av, alpha_av*tau_av, tau_at, alpha_at*tau_at);
fprintf('   Forma Polo-Zero    : D(s) = %.6f * [(s + %.2f)(s + %.2f)] / [(s + %.2f)(s + %.6f)]\n', ...
    Kc_b, z_av, z_at, p_av, p_at);
fprintf('   Forma Polinomial   : (%.6f*s^2 + %.6f*s + %.6f) / (%.6f*s^2 + %.6f*s + %.6f)\n\n', ...
    num_Db(1), num_Db(2), num_Db(3), den_Db(1), den_Db(2), den_Db(3));

fprintf('4. Verificação das Especificações (Item b):\n');
fprintf('   Kv obtido = %.4f s^-1 (Especificado: %g s^-1) -> %s\n', Kv_b, Kv_des, status_Kv_b);
fprintf('   MG obtida = %.4f dB em w_cf = %.4f rad/s (Especificado: >= %g dB) -> %s\n', ...
    MG_b_dB, wcf_b, MG_des_dB, status_MG_b);
fprintf('   MF obtida = %.4f deg em w_cg = %.4f rad/s (Especificado: >= %g deg) -> %s\n', ...
    MF_b, wcg_b, MF_des, status_MF_b);
fprintf('   Polos de Malha Fechada T_b(s):\n');
disp(polos_Tb);
fprintf('   Desempenho Transitório T_b(s): tr = %.4f s | Mp = %.4f %% | ts (2%%) = %.4f s\n\n', ...
    info_b.RiseTime, info_b.Overshoot, info_b.SettlingTime);

% -------------------------------------------------------------------------
% Gráfico (b): Bode Compensado (mesmo padrão visual do item a.2)
% -------------------------------------------------------------------------
fig_b = figure('Name', 'Item (b) - Bode Compensado', 'Position', [180 120 860 680]);
w_bode_b = logspace(-2.5, 2.4, 3000);
[mag_b_raw, phase_b_raw] = bode(L_b, w_bode_b);
[mag_g0_raw, phase_g0_raw] = bode(G, w_bode_b);
mag_b_dB = 20 * log10(squeeze(mag_b_raw));
mag_g0_dB = 20 * log10(squeeze(mag_g0_raw));
phase_b_deg = unwrap(squeeze(phase_b_raw) * pi/180) * 180/pi;
phase_g0_deg = unwrap(squeeze(phase_g0_raw) * pi/180) * 180/pi;
if phase_b_deg(1) > 0, phase_b_deg = phase_b_deg - 360; end
if phase_g0_deg(1) > 0, phase_g0_deg = phase_g0_deg - 360; end

% Subgráfico Superior: Magnitude (dB)
ax_b_mag = subplot(2, 1, 1);
semilogx(w_bode_b, mag_b_dB, '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2, 'DisplayName', 'Compensado D(s)G(s)'); hold on;
semilogx(w_bode_b, mag_g0_dB, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3, 'DisplayName', 'Não Compensado G(s)');
yline(0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.1, 'HandleVisibility', 'off');

plot(wcg_b, 0, 'o', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 8, 'HandleVisibility', 'off');
text(wcg_b * 1.25, 10, sprintf('\\omega_{cg} = %.2f rad/s (0 dB)', wcg_b), ...
    'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 10);

plot([wcf_b wcf_b], [mag_wcf_b_dB 0], '-', 'Color', [0.55 0 0.65], 'LineWidth', 2.2, 'HandleVisibility', 'off');
plot(wcf_b, mag_wcf_b_dB, 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8, 'HandleVisibility', 'off');
text(wcf_b * 1.20, mag_wcf_b_dB / 2, sprintf('MG = %.2f dB', MG_b_dB), ...
    'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 10);

xlim([10^-2.5 250]); ylim([-100 70]);
xticks([1e-2 1e-1 1e0 1e1 1e2]); yticks(-100:25:75);
ylabel('Magnitude (dB)', 'FontSize', 11);
title(sprintf('Diagrama de Bode Compensado D(s)G(s) | Kv = %.1f s^{-1} | MG = %.2f dB | MF = %.2f^{\\circ}', ...
    Kv_b, MG_b_dB, MF_b), 'FontSize', 12.5, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 9.5);
grid on; grid minor; box on;

% Subgráfico Inferior: Fase (graus)
ax_b_pha = subplot(2, 1, 2);
semilogx(w_bode_b, phase_b_deg, '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2); hold on;
semilogx(w_bode_b, phase_g0_deg, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3);
yline(-180, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);

plot([wcg_b wcg_b], [-180 phase_wcg_b], '-', 'Color', [0.8 0 0], 'LineWidth', 2.2);
plot(wcg_b, phase_wcg_b, 'o', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 8);
text(wcg_b * 1.25, (-180 + phase_wcg_b) / 2, sprintf('MF = %.2f^{\\circ}', MF_b), ...
    'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 10);

plot(wcf_b, -180, 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8);
text(wcf_b * 1.20, -170, sprintf('\\omega_{cf} = %.2f rad/s (-180^{\\circ})', wcf_b), ...
    'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 10);

xlim([10^-2.5 250]); ylim([-278 -80]);
xticks([1e-2 1e-1 1e0 1e1 1e2]); yticks(-270:45:-90);
xlabel('Frequência (rad/s)', 'FontSize', 11);
ylabel('Fase (graus)', 'FontSize', 11);
grid on; grid minor; box on;
linkaxes([ax_b_mag, ax_b_pha], 'x');
aplicar_light_mode_e_exportar(fig_b, 'item_b_Bode_Compensado.png');


%% ========================================================================
%  ITEM (c): CONTROLADOR POR ALOCAÇÃO DE POLOS (FUNÇÃO DE TRANSFERÊNCIA / LGR)
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (c): CONTROLADOR POR ALOCAÇÃO DE POLOS (JUSTIFICADO PELO LGR)\n');
fprintf('*************************************************************************\n');

% Polos desejados escolhidos para obter resposta transitória semelhante ao item (b)
% e garantir que o controlador C_c(s) = (b1*s + b0)/(s + a0) seja de fase mínima.
s1_des = -0.90 + 1j * 0.50;
s2_des = conj(s1_des);
s3_des = -7.00;
wn_c   = abs(s1_des);
zeta_c = -real(s1_des) / wn_c;

% Polinômio parcial P3(s) = (s - s1)*(s - s2)*(s - s3) = s^3 + c2*s^2 + c1*s + c0
P3 = real(poly([s1_des, s2_des, s3_des]));
c2 = P3(2); c1 = P3(3); c0 = P3(4);

% Equação Característica: 1 + C_c(s)*G(s) = 0
% s*(s+8)*(s+20)*(s+a0) + 160*(b1*s + b0) = 0
% s^4 + (28 + a0)*s^3 + (160 + 28*a0)*s^2 + (160*a0 + 160*b1)*s + 160*b0 = 0
s4_des = (160 + 28*c2 - 784 - c1) / (28 - c2);
a0_c = c2 - 28 - s4_des;

P4_des = real(conv(P3, [1 -s4_des])); % s^4 + alpha3*s^3 + alpha2*s^2 + alpha1*s + alpha0
alpha3_c = P4_des(2); alpha2_c = P4_des(3); alpha1_c = P4_des(4); alpha0_c = P4_des(5);

b1_c = (alpha1_c - 160 * a0_c) / 160;
b0_c = alpha0_c / 160;
zero_Cc = -b0_c / b1_c;
polo_Cc = -a0_c;

C_c = tf([b1_c b0_c], [1 a0_c]);
L_c = series(C_c, G);
T_c = feedback(L_c, 1);

Kv_c = dcgain(s * L_c);
[MG_c_abs, MF_c, wcf_c, wcg_c] = margin(L_c);
MG_c_dB = 20 * log10(MG_c_abs);
polos_Tc = pole(T_c);
info_c = stepinfo(T_c);

fprintf('1. Especificação dos Polos Desejados:\n');
fprintf('   zeta = %.4f | wn = %.4f rad/s\n', zeta_c, wn_c);
fprintf('   Polos dominantes : s1,2 = %.4f +- j%.4f\n', real(s1_des), imag(s1_des));
fprintf('   Polos auxiliares : s3 = %.4f | s4 (resultante) = %.4f\n\n', s3_des, s4_des);

fprintf('2. Coeficientes Intermediários dos Polinômios:\n');
fprintf('   Polinômio parcial P3(s): c2 = %.4f | c1 = %.4f | c0 = %.4f\n', c2, c1, c0);
fprintf('   Polinômio desejado P4(s): alpha3 = %.4f | alpha2 = %.4f | alpha1 = %.4f | alpha0 = %.4f\n', ...
    alpha3_c, alpha2_c, alpha1_c, alpha0_c);
fprintf('   Parâmetros de C_c(s): a0 = %.6f | b1 = %.6f | b0 = %.6f\n\n', a0_c, b1_c, b0_c);

fprintf('3. Função de Transferência C_c(s) = (%.6f*s + %.6f) / (s + %.6f)\n', b1_c, b0_c, a0_c);
fprintf('   Zero do controlador: s = %.4f (Fase Mínima: %s)\n', zero_Cc, verificar_status(zero_Cc < 0));
fprintf('   Polo do controlador: s = %.4f (Estável: %s)\n\n', polo_Cc, verificar_status(polo_Cc < 0));

fprintf('4. Verificação de Desempenho (Item c):\n');
fprintf('   Polos de Malha Fechada T_c(s):\n');
disp(polos_Tc);
fprintf('   Desempenho Transitório T_c(s): tr = %.4f s | Mp = %.4f %% | ts (2%%) = %.4f s\n', ...
    info_c.RiseTime, info_c.Overshoot, info_c.SettlingTime);
fprintf('   Margens e Erro: Kv = %.4f s^-1 | MG = %.4f dB | MF = %.4f deg\n\n', ...
    Kv_c, MG_c_dB, MF_c);

% -------------------------------------------------------------------------
% Gráfico (c): LGR do Sistema Compensado por Alocação de Polos
% -------------------------------------------------------------------------
fig_c = figure('Name', 'Item (c) - LGR Alocação de Polos', 'Position', [220 120 840 660]);
k_lgr_c = unique([0, logspace(-3, 2.5, 4000), 1]);
r_lgr_c = rlocus(L_c, k_lgr_c);
polos_Lc = pole(L_c);
zeros_Lc = zero(L_c);

hold on;
xline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');

plot(real(r_lgr_c(1,:)), imag(r_lgr_c(1,:)), '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2, 'DisplayName', 'Ramos LGR de C_c(s)G(s)');
for idx_r = 2:size(r_lgr_c, 1)
    plot(real(r_lgr_c(idx_r,:)), imag(r_lgr_c(idx_r,:)), '-', 'Color', [0 0.28 0.72], 'LineWidth', 2.2, 'HandleVisibility', 'off');
end

plot(real(polos_Lc), imag(polos_Lc), 'rx', 'MarkerSize', 9, 'LineWidth', 2.2, 'DisplayName', 'Polos de Malha Aberta');
plot(real(zeros_Lc), imag(zeros_Lc), 'bo', 'MarkerSize', 9, 'LineWidth', 2.2, 'DisplayName', 'Zero de C_c(s)');
plot(real(polos_Tc), imag(polos_Tc), 's', 'Color', [0.55 0 0.65], 'MarkerFaceColor', [0.55 0 0.65], 'MarkerSize', 8.5, 'DisplayName', 'Polos Alocados em MF (K=1)');

text(real(s1_des) + 0.4, imag(s1_des) + 1.1, sprintf('s_{1,2} = %.2f \\pm j%.2f', real(s1_des), imag(s1_des)), ...
    'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 9.5);
text(s3_des, 1.3, sprintf('s_3 = %.2f', s3_des), ...
    'HorizontalAlignment', 'center', 'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 9.5);
text(s4_des, 1.3, sprintf('s_4 = %.2f', s4_des), ...
    'HorizontalAlignment', 'center', 'Color', [0.55 0 0.65], 'FontWeight', 'bold', 'FontSize', 9.5);
text(zero_Cc, -1.4, sprintf('Zero (%.2f)', zero_Cc), ...
    'HorizontalAlignment', 'center', 'Color', [0 0.28 0.72], 'FontWeight', 'bold', 'FontSize', 9.5);
text(polo_Cc, -1.4, sprintf('Polo C_c (%.2f)', polo_Cc), ...
    'HorizontalAlignment', 'center', 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 9.5);

xlim([-27 4]); ylim([-14 14]);
xticks(-25:5:5); yticks(-12:4:12);
title('Lugar das Raízes de C_c(s)G(s) com Polos Alocados em Malha Fechada', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Eixo Real (s^{-1})', 'FontSize', 11);
ylabel('Eixo Imaginário (rad/s)', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9.5);
grid on; grid minor; box off;
aplicar_light_mode_e_exportar(fig_c, 'item_c_LGR_Alocacao_Polos.png');


%% ========================================================================
%  ITEM (d): CONTROLADOR DO TIPO PID (PID POR ALOCAÇÃO DE POLOS)
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (d): PROJETO DE CONTROLADOR PID POR ALOCAÇÃO DE POLOS\n');
fprintf('*************************************************************************\n');

% Conforme slide "PID por Alocação de Polos": 1 + Gc(s)*G(s) = prod(s - pi) = 0
% Gc(s) = Kp*(1 + 1/(Ti*s) + Td*s) = Kp + Ki/s + Kd*s = (Kd*s^2 + Kp*s + Ki) / s
% Eq. Característica: s^2*(s+8)*(s+20) + 160*(Kd*s^2 + Kp*s + Ki) = 0
% s^4 + 28*s^3 + (160 + 160*Kd)*s^2 + 160*Kp*s + 160*Ki = 0
% Como a planta é de 3ª ordem (MF de 4ª ordem) e o PID possui 3 parâmetros (Kp, Ki, Kd),
% a soma dos 4 polos é invariante (-28), vinculando o 4º polo aos 3 escolhidos.
s1_pid = -0.55;
s2_pid = -0.85;
s3_pid = -8.00;
s4_pid = -28.0 - (s1_pid + s2_pid + s3_pid); % Imposto pelo coeficiente de s^3 (-28)

P4_pid = real(poly([s1_pid, s2_pid, s3_pid, s4_pid]));
alpha3_pid = P4_pid(2); % Exatamente 28.0000
alpha2_pid = P4_pid(3);
alpha1_pid = P4_pid(4);
alpha0_pid = P4_pid(5);

Kd = (alpha2_pid - 160) / 160;
Kp = alpha1_pid / 160;
Ki = alpha0_pid / 160;
Ti_pid = Kp / Ki;
Td_pid = Kd / Kp;

zeros_PID = roots([Kd Kp Ki]);
C_pid = tf([Kd Kp Ki], [1 0]);
L_pid = series(C_pid, G);
T_d = feedback(L_pid, 1);

[MG_d_abs, MF_d, wcf_d, wcg_d] = margin(L_pid);
MG_d_dB = 20 * log10(MG_d_abs);
polos_Td = pole(T_d);
info_d = stepinfo(T_d);

fprintf('1. Polos Desejados e Coeficientes da Equação Característica (Grau 4):\n');
fprintf('   Polos de MF escolhidos : s1 = %.4f, s2 = %.4f, s3 = %.4f\n', s1_pid, s2_pid, s3_pid);
fprintf('   4o Polo (traço = -28)  : s4 = %.4f\n', s4_pid);
fprintf('   Coeficientes desejados : alpha3 = %.4f | alpha2 = %.4f | alpha1 = %.4f | alpha0 = %.4f\n\n', ...
    alpha3_pid, alpha2_pid, alpha1_pid, alpha0_pid);

fprintf('2. Parâmetros Calculados do Controlador PID Gc(s) = Kp*(1 + 1/(Ti*s) + Td*s):\n');
fprintf('   Kp = %.6f\n', Kp);
fprintf('   Ki = %.6f  (Ti = Kp/Ki = %.6f s)\n', Ki, Ti_pid);
fprintf('   Kd = %.6f  (Td = Kd/Kp = %.6f s)\n', Kd, Td_pid);
fprintf('   Zeros introduzidos pelo PID: z1 = %.4f | z2 = %.4f\n\n', zeros_PID(1), zeros_PID(2));

fprintf('3. Verificação de Desempenho (Item d):\n');
fprintf('   Polos de Malha Fechada T_d(s):\n');
disp(polos_Td);
fprintf('   Desempenho Transitório T_d(s): tr = %.4f s | Mp = %.4f %% | ts (2%%) = %.4f s\n', ...
    info_d.RiseTime, info_d.Overshoot, info_d.SettlingTime);
fprintf('   Tipo do Sistema: Tipo 2 -> Kv = Inf (erro nulo à rampa) | MG = %.4f dB | MF = %.4f deg\n\n', ...
    MG_d_dB, MF_d);


%% ========================================================================
%  ITEM (e): CONTROLADOR POR REALIMENTAÇÃO DE ESTADO (AULA 17)
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (e): PROJETO POR REALIMENTAÇÃO DE ESTADO (u(t) = -K*x(t) + J*r(t))\n');
fprintf('*************************************************************************\n');

% Forma Canônica Controlável (Variável de Fase - Aula 17, Slide 10 e 14)
A = [0    1    0;
     0    0    1;
     0 -160  -28];
B = [0; 0; 1];
C = [160 0 0];
D = 0;

% Matriz de Controlabilidade S = [B  A*B  A^2*B] (Aula 17, Slide 10)
S = ctrb(A, B);
posto_S = rank(S);
det_S = det(S);

% Polos desejados de malha fechada (mesmos polos dominantes do item (c))
polos_des_e = [s1_des, s2_des, -12.0];
P3_e = real(poly(polos_des_e));
alpha2_e = P3_e(2);
alpha1_e = P3_e(3);
alpha0_e = P3_e(4);

% Cálculo analítico do vetor de ganhos K = [k1 k2 k3] (Aula 17, Slide 18)
% EC: det(sI - (A - B*K)) = s^3 + (28 + k3)*s^2 + (160 + k2)*s + k1
k1_sf = alpha0_e;
k2_sf = alpha1_e - 160;
k3_sf = alpha2_e - 28;
K_sf = [k1_sf, k2_sf, k3_sf];

% Conferência com comando MATLAB place(A, B, P) (Aula 17, Slide 19)
K_place = place(A, B, polos_des_e);

% Ganho J tal que o ganho CC de r(t) para y(t) garanta erro nulo ao degrau:
% T(s->0) = C*(-A + B*K)^-1 * B * J = 1  (Aula 17, Slide 19)
A_bar = A - B * K_sf;
MJ = C * ((-A + B * K_sf) \ B);
J = 1 / MJ;

sys_e = ss(A_bar, B * J, C, D);
T_e = tf(sys_e);
polos_Te = eig(A_bar);
info_e = stepinfo(sys_e);

fprintf('1. Representação em Espaço de Estados e Controlabilidade:\n');
fprintf('   det(S) = %.4f | Posto(S) = %d (Controlável: %s)\n\n', ...
    det_S, posto_S, verificar_status(posto_S == size(A, 1)));

fprintf('2. Equação Característica Desejada e Ganhos de Realimentação:\n');
fprintf('   Polos desejados: %.4f +- j%.4f e %.4f\n', real(s1_des), imag(s1_des), polos_des_e(3));
fprintf('   Coeficientes   : alpha2 = %.4f | alpha1 = %.4f | alpha0 = %.4f\n', alpha2_e, alpha1_e, alpha0_e);
fprintf('   Vetor K (Analítico): [k1, k2, k3] = [%.6f, %.6f, %.6f]\n', K_sf(1), K_sf(2), K_sf(3));
fprintf('   Vetor K (place)    : [k1, k2, k3] = [%.6f, %.6f, %.6f]\n', K_place(1), K_place(2), K_place(3));
fprintf('   Ganho de referência J (para erro nulo ao degrau): %.6f\n\n', J);

fprintf('3. Verificação de Desempenho e Observações Práticas (Item e):\n');
fprintf('   Autovalores de A_barra = (A - B*K):\n');
disp(polos_Te);
fprintf('   Desempenho Transitório: tr = %.4f s | Mp = %.4f %% | ts (2%%) = %.4f s\n', ...
    info_e.RiseTime, info_e.Overshoot, info_e.SettlingTime);
fprintf('   Ganho CC de malha fechada T_e(0) = %.4f (erro ao degrau ess = %.4f)\n', ...
    dcgain(sys_e), abs(1 - dcgain(sys_e)));
fprintf('   Nota Técnica: k2 (%.2f) e k3 (%.2f) são negativos porque os polos de malha aberta\n', k2_sf, k3_sf);
fprintf('   (-8 e -20) são mais rápidos que os polos dominantes desejados (-0.90 +- j0.50), exigindo\n');
fprintf('   redução do amortecimento interno da planta. Como x2 e x3 correspondem às derivadas da\n');
fprintf('   saída (y_dot/160 e y_ddot/160), a implementação física requer um observador de estados.\n\n');


%% ========================================================================
%  ITEM (f): RESPOSTA AO DEGRAU, TABELA COMPARATIVA E ANÁLISE DINÂMICA
%% ========================================================================
fprintf('*************************************************************************\n');
fprintf('ITEM (f): COMPARAÇÃO DAS RESPOSTAS AO DEGRAU UNITÁRIO E TABELA RESUMO\n');
fprintf('*************************************************************************\n');

ess_b = abs(1 - dcgain(T_b));
ess_c = abs(1 - dcgain(T_c));
ess_d = abs(1 - dcgain(T_d));
ess_e = abs(1 - dcgain(sys_e));

fprintf('+-----------------------------------+------------+---------------+-----------------+-------------+------------+\n');
fprintf('| Controlador                       | tr (s)     | Sobrepasso(%%) | ts - 2%% (s)     | ess (Degrau)| Kv (s^-1)  |\n');
fprintf('+-----------------------------------+------------+---------------+-----------------+-------------+------------+\n');
fprintf('| (b) Freq. (Avanço-Atraso)         | %10.4f | %13.4f | %15.4f | %11.4f | %10.4f |\n', ...
    info_b.RiseTime, info_b.Overshoot, info_b.SettlingTime, ess_b, Kv_b);
fprintf('| (c) Alocação de Polos (FT/LGR)    | %10.4f | %13.4f | %15.4f | %11.4f | %10.4f |\n', ...
    info_c.RiseTime, info_c.Overshoot, info_c.SettlingTime, ess_c, Kv_c);
fprintf('| (d) Controlador PID               | %10.4f | %13.4f | %15.4f | %11.4f |        Inf |\n', ...
    info_d.RiseTime, info_d.Overshoot, info_d.SettlingTime, ess_d);
fprintf('| (e) Realimentação de Estado       | %10.4f | %13.4f | %15.4f | %11.4f |   N/A (T0) |\n', ...
    info_e.RiseTime, info_e.Overshoot, info_e.SettlingTime, ess_e);
fprintf('+-----------------------------------+------------+---------------+-----------------+-------------+------------+\n\n');

ts_vetor = [info_b.SettlingTime, info_c.SettlingTime, info_d.SettlingTime, info_e.SettlingTime];
mp_vetor = [info_b.Overshoot,    info_c.Overshoot,    info_d.Overshoot,    info_e.Overshoot];
nomes_ctrl = {'(b) Frequência', '(c) Alocação de Polos', '(d) PID', '(e) Realimentação de Estado'};
[~, idx_menor_ts] = min(ts_vetor);
[~, idx_menor_mp] = min(mp_vetor);

fprintf('ANÁLISE COMPARATIVA (BASEADA NOS VALORES CALCULADOS ACIMA):\n');
fprintf('1. Comparação de Desempenho Transitório e Estático:\n');
fprintf('   - Tempo de subida (tr): varia de %.2f s no PID (d) até %.2f s na Realimentação de Estado (e),\n', ...
    info_d.RiseTime, info_e.RiseTime);
fprintf('     com valores intermediários de %.2f s em (c) e %.2f s em (b).\n', ...
    info_c.RiseTime, info_b.RiseTime);
fprintf('   - Sobrepasso (Mp): o menor sobrepasso ocorre em %s (%.2f%%). Os controladores (b) (%.2f%%)\n', ...
    nomes_ctrl{idx_menor_mp}, mp_vetor(idx_menor_mp), info_b.Overshoot);
fprintf('     e (d) (%.2f%%) exibem sobrepasso moderado devido aos zeros próximos da origem introduzidos\n', ...
    info_d.Overshoot);
fprintf('     pela rede de atraso (z_at = -%.2f) e pelo integrador do PID (z1 = %.3f), respectivamente.\n', ...
    z_at, zeros_PID(1));
fprintf('   - Tempo de assentamento (ts): o mais rápido é %s (%.2f s). Já o controlador (b) apresenta\n', ...
    nomes_ctrl{idx_menor_ts}, ts_vetor(idx_menor_ts));
fprintf('     cauda lenta (ts = %.2f s) provocada pelo dipolo polo-zero próximo à origem (p_at = -%.5f,\n', ...
    info_b.SettlingTime, p_at);
fprintf('     z_at = -%.2f) necessário para elevar Kv a %.1f s^-1 mantendo MG = %.2f dB e MF = %.2f deg.\n', ...
    z_at, Kv_b, MG_b_dB, MF_b);
fprintf('   - Erro em regime permanente: todos anulam o erro ao degrau (ess = 0). Para entrada em rampa,\n');
fprintf('     (b) atende exatamente Kv = %.2f s^-1, (c) reduz Kv para %.2f s^-1 e o PID (d) eleva o sistema\n', ...
    Kv_b, Kv_c);
fprintf('     para Tipo 2 (Kv = Inf, erro nulo também à rampa).\n\n');

fprintf('2. Escolha do Melhor Controlador:\n');
fprintf('   - Para atender rigorosamente às especificações do enunciado (Kv = %g s^-1, MG >= %g dB e MF >= %g deg),\n', ...
    Kv_des, MG_des_dB, MF_des);
fprintf('     o Controlador no Domínio da Frequência (item b) é o mais indicado, pois garante simultaneamente\n');
fprintf('     precisão de velocidade e margens de estabilidade robustas (MG = %.2f dB, MF = %.2f deg).\n', ...
    MG_b_dB, MF_b);
fprintf('   - Se o critério principal for a qualidade transitória pura (menor sobrepasso = %.2f%% e assentamento\n', ...
    info_e.Overshoot);
fprintf('     rápido = %.2f s sem zeros finitos), a Realimentação de Estado (item e) oferece o melhor comportamento\n', ...
    info_e.SettlingTime);
fprintf('     dinâmico, sendo recomendada na prática com adição de ação integral e observador de estados.\n\n');

% -------------------------------------------------------------------------
% Gráfico (f): Resposta ao Degrau Unitário Comparativa
% -------------------------------------------------------------------------
t_sim = 0:0.02:22;
[y_b, t_b] = step(T_b, t_sim);
[y_c, t_c] = step(T_c, t_sim);
[y_d, t_d] = step(T_d, t_sim);
[y_e, t_e] = step(sys_e, t_sim);

fig_f = figure('Name', 'Item (f) - Comparação Resposta ao Degrau', 'Position', [260 120 900 580]);
hold on;
% Faixa de tolerância de +-2% e referência unitária
fill([0 22 22 0], [0.98 0.98 1.02 1.02], [0.92 0.94 0.98], 'EdgeColor', 'none', 'DisplayName', 'Faixa de Acomodação (\pm2%)');
yline(1.0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2, 'DisplayName', 'Referência (Degrau Unitário)');

plot(t_b, y_b, '-',  'Color', [0 0.28 0.72],  'LineWidth', 2.2, 'DisplayName', '(b) Domínio da Frequência (Avanço-Atraso)');
plot(t_c, y_c, '--', 'Color', [0.82 0.1 0.1], 'LineWidth', 2.2, 'DisplayName', '(c) Alocação de Polos (FT/LGR)');
plot(t_d, y_d, '-.', 'Color', [0 0.55 0.22],  'LineWidth', 2.2, 'DisplayName', '(d) Controlador PID');
plot(t_e, y_e, ':',  'Color', [0.55 0 0.65],  'LineWidth', 2.5, 'DisplayName', '(e) Realimentação de Estado');

legend('Location', 'southeast', 'FontSize', 9.5);
title('Resposta ao Degrau Unitário em Malha Fechada (Itens b a e)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Amplitude y(t)', 'FontSize', 11);
xlim([0 22]); ylim([0 1.18]);
xticks(0:2:22); yticks(0:0.2:1.2);
grid on; grid minor; box off;
aplicar_light_mode_e_exportar(fig_f, 'item_f_Comparacao_Degrau.png');

%% ========================================================================
%  FUNÇÕES AUXILIARES
%% ========================================================================
function txt = verificar_status(condicao)
    if condicao
        txt = 'ATENDIDO';
    else
        txt = 'NÃO ATENDIDO';
    end
end

function aplicar_light_mode_e_exportar(fig, nome_arquivo)
    % Compatível com MATLAB R2024b/R2025a/R2026a+ (tema escuro nativo)
    try
        theme(fig, 'light');
    catch
        % Ignora caso a versão do MATLAB seja anterior ao suporte de theme()
    end

    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');

    eixos = findall(fig, 'Type', 'axes');
    for k = 1:length(eixos)
        set(eixos(k), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
            'LineWidth', 0.9, 'TickDir', 'in', 'FontSize', 10.5, ...
            'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, ...
            'MinorGridColor', [0.90 0.90 0.90], 'MinorGridAlpha', 1, ...
            'MinorGridLineStyle', ':');
        if ~isempty(eixos(k).Title)
            set(eixos(k).Title, 'Color', 'k');
        end
        if ~isempty(eixos(k).XLabel)
            set(eixos(k).XLabel, 'Color', 'k');
        end
        if ~isempty(eixos(k).YLabel)
            set(eixos(k).YLabel, 'Color', 'k');
        end
    end

    legendas = findall(fig, 'Type', 'legend');
    for k = 1:length(legendas)
        set(legendas(k), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.75 0.75 0.75]);
    end

    drawnow;
    exportgraphics(fig, nome_arquivo, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('>> Gráfico exportado em Light Mode: %s\n\n', nome_arquivo);
end