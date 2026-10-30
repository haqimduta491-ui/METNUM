%% ============================================================
%  RESPONSI 1 PRAKTIKUM METODE NUMERIK
%  Kasus 3 : Persamaan Nonlinear (Bola di dalam Tabung)
%% ============================================================
clc; clear; close all;

%% ---------- 1. DATA DAN PERSAMAAN NONLINEAR ----------
% Tabung: t = 2r -> V_tabung = pi r^2 (2r) = 2 pi r^3
% V_bola : V_tabung = (4/3) pi r^3 : 2 pi r^3 = 2 : 3
V_bola   = 300;
V_tabung = 1.5 * V_bola;
V_luar   = V_tabung - V_bola;
f = @(r) (2/3)*pi*r.^3 - V_luar;

%% ---------- 2. METODE GRAFIK ----------
r = linspace(0, 8, 401);
figure(1);
plot(r, f(r), 'b-', 'LineWidth', 2); hold on;
plot([0 8], [0 0], 'k-');
plot(4, f(4), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(5, f(5), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
grid on;
xlabel('r (cm)'); ylabel('f(r)');
title('Grafik f(r) = (2/3)\pi r^3 - 150');
legend('f(r)', 'sumbu r', 'a = 4', 'b = 5', 'Location', 'northwest');
hold off;
fprintf('f(4) = %.4f\n', f(4));
fprintf('f(5) = %.4f\n', f(5));

%% ---------- 3. METODE BISEKSI ----------
a    = 4;
b    = 5;
tol  = 1e-6;
maks = 100;
if f(a)*f(b) >= 0
    error('Tidak ada akar pada selang [a, b]');
end
fprintf('\nIterasi Biseksi:\n');
fprintf('Iter |      a           b           c          f(c)      | Galat\n');
for iter = 1:maks
    c  = (a + b) / 2;
    fc = f(c);
    galat = (b - a) / 2;
    fprintf('%4d | %10.6f %10.6f %10.6f %12.6f | %.3e\n', iter, a, b, c, fc, galat);
    if f(a)*fc < 0
        b = c;
    else
        a = c;
    end
    if galat < tol
        break;
    end
end

%% ---------- 4. HASIL AKHIR ----------
r_bola = c;
t_tabung = 2 * r_bola;
V_bola_cek   = (4/3)*pi*r_bola^3;
V_tabung_cek = pi*r_bola^2*t_tabung;
V_luar_cek   = V_tabung_cek - V_bola_cek;
r_analitik   = (3*V_luar/(2*pi))^(1/3);
fprintf('\nHasil akhir:\n');
fprintf('Jari-jari bola dan tabung : %.6f cm\n', r_bola);
fprintf('Tinggi tabung (2r)        : %.6f cm\n', t_tabung);
fprintf('Volume bola               : %.4f cm^3\n', V_bola_cek);
fprintf('Volume tabung             : %.4f cm^3\n', V_tabung_cek);
fprintf('Volume tabung di luar bola: %.4f cm^3\n', V_luar_cek);
fprintf('Jumlah iterasi            : %d\n', iter);
fprintf('r analitik (pembanding)   : %.6f cm\n', r_analitik);
fprintf('Selisih dengan analitik   : %.3e cm\n', abs(r_bola - r_analitik));
