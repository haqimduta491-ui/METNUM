% =========================================================
% Perbandingan Metode Euler, Heun (RK2), dan RK4
% PDB : dy/dx = x*sqrt(y)
% Asumsi nilai awal : x0 = 0, y0 = 1
% Target            : y(1.2) dengan h = 0.1
% Toleransi error   : eps < 0.0001
% =========================================================
clear; clc;

% ---------- Definisi masalah ----------
f      = @(x, y) x .* sqrt(y);          % f(x,y)
yexact = @(x) (x.^2 / 4 + 1).^2;        % solusi eksak
x0 = 0;  y0 = 1;
h  = 0.1;
xT = 1.2;
tol = 0.0001;
N  = round((xT - x0) / h);              % jumlah langkah (=12)

% ---------- Penyimpanan hasil ----------
x  = x0 + (0:N) * h;
yE = zeros(1, N+1);   % Euler
yH = zeros(1, N+1);   % Heun
yR = zeros(1, N+1);   % RK4
yE(1) = y0;  yH(1) = y0;  yR(1) = y0;

% ---------- Iterasi ----------
for i = 1:N
    xi = x(i);

    % 1. Euler
    yE(i+1) = yE(i) + h * f(xi, yE(i));

    % 2. Heun (Runge-Kutta orde-2)
    k1 = f(xi, yH(i));
    k2 = f(xi + h, yH(i) + h * k1);
    yH(i+1) = yH(i) + (h/2) * (k1 + k2);

    % 3. Runge-Kutta orde-4
    k1 = f(xi, yR(i));
    k2 = f(xi + 0.5*h, yR(i) + 0.5*h*k1);
    k3 = f(xi + 0.5*h, yR(i) + 0.5*h*k2);
    k4 = f(xi + h,     yR(i) + h*k3);
    yR(i+1) = yR(i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

yEx = yexact(x);

% ---------- Tabel semua langkah ----------
printf('\nTABEL ITERASI (h = %.1f)\n', h);
printf('%-5s %-6s %-12s %-12s %-12s %-12s\n', ...
       'i', 'x', 'Euler', 'Heun', 'RK4', 'Eksak');
for i = 1:N+1
    printf('%-5d %-6.1f %-12.6f %-12.6f %-12.6f %-12.6f\n', ...
           i-1, x(i), yE(i), yH(i), yR(i), yEx(i));
end

% ---------- Perbandingan di x = 1.2 ----------
hasil = [yE(end), yH(end), yR(end)];
nama  = {'Euler', 'Heun (RK2)', 'RK4'};
err   = abs(hasil - yEx(end));

printf('\nPERBANDINGAN PADA x = %.1f\n', xT);
printf('Solusi eksak y(%.1f) = %.8f\n\n', xT, yEx(end));
printf('%-12s %-14s %-14s %-10s\n', 'Metode', 'y(1.2)', 'Error mutlak', 'Error<tol?');
for j = 1:3
    if err(j) < tol
        status = 'YA';
    else
        status = 'TIDAK';
    end
    printf('%-12s %-14.8f %-14.8f %-10s\n', nama{j}, hasil(j), err(j), status);
end

% ---------- Grafik solusi ----------
xf = linspace(x0, xT, 200);
figure(1);
plot(xf, yexact(xf), 'k-', 'LineWidth', 2); hold on;
plot(x, yE, 'ro--', x, yH, 'gs--', x, yR, 'b^--');
hold off; grid on;
xlabel('x'); ylabel('y');
title('dy/dx = x\surd y, y(0)=1, h=0.1');
legend('Eksak', 'Euler', 'Heun', 'RK4', 'Location', 'northwest');

% ---------- Grafik error ----------
figure(2);
semilogy(x(2:end), abs(yE(2:end)-yEx(2:end)), 'ro-', ...
         x(2:end), abs(yH(2:end)-yEx(2:end)), 'gs-', ...
         x(2:end), abs(yR(2:end)-yEx(2:end)) + eps, 'b^-');
hold on; semilogy([x(2) x(end)], [tol tol], 'k--'); hold off;
grid on; xlabel('x'); ylabel('Error mutlak (skala log)');
title('Error terhadap solusi eksak');
legend('Euler', 'Heun', 'RK4', 'Toleransi 10^{-4}', 'Location', 'southeast');
