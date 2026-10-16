
clc; clear all; close all;



f  = @(x) log10(x+2) + 2*x - 3;
g1 = @(x) log10(x+2);
g2 = @(x) 3 - 2*x;

% Rentang sumbu x
xmin = -1.9;
xmax = 3;

% Judul dan keterangan
judul = 'Persamaan 3';
ket1  = 'y = log(x+2)';
ket2  = 'y = 3 - 2x';
% =====================================================================

x = linspace(xmin, xmax, 1000);

% ---------------------- GRAFIK TUNGGAL ----------------------
figure(1);
plot(x, f(x), 'b-', 'LineWidth', 2); hold on;
plot(x, zeros(size(x)), 'k--', 'LineWidth', 1);
grid on;
xlabel('x'); ylabel('f(x)');
title([judul ' - Grafik Tunggal']);
legend('f(x)', 'y = 0', 'Location', 'best');

% ---------------------- GRAFIK GANDA ------------------------
figure(2);
plot(x, g1(x), 'b-', 'LineWidth', 2); hold on;
plot(x, g2(x), 'r-', 'LineWidth', 2);
grid on;
xlabel('x'); ylabel('y');
title([judul ' - Grafik Ganda']);
legend(ket1, ket2, 'Location', 'best');
