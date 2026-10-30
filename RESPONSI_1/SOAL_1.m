%% ============================================================
%  RESPONSI 1 PRAKTIKUM METODE NUMERIK
%  Kasus 1 : Galat & Error Pemotongan (Gelombang Air Laut)
%  h = h0 * ( cos(t*v/lambda) + e^x )
%% ============================================================
clc; clear; close all;

%% ---------- 1. INPUT PARAMETER ----------
lambda = 12;                 % panjang gelombang
v      = 36;                 % kecepatan
x      = 10;                 % nilai x pada e^x
h0     = 5;                  % amplitudo
t      = linspace(0, 5, 501);% t dari 0 s.d. 5 (501 titik)
n_list = [23 24 25];         % batas iterasi deret
tol    = 1;                  % toleransi error (%)

theta = t * v / lambda;      % argumen cos: tv/lambda

%% ---------- 2. HASIL EKSAK ----------
h_eksak = h0 * (cos(theta) + exp(x));

%% ---------- 3. HASIL PENDEKATAN (DERET MACLAURIN) ----------
h_approx = zeros(length(n_list), length(t));   % baris = tiap n

for idx = 1:length(n_list)
    n = n_list(idx);

    % --- pendekatan cos(theta) = sum_{k=0}^{n} (-1)^k theta^(2k)/(2k)! ---
    cos_n = zeros(size(theta));
    for k = 0:n
        cos_n = cos_n + ((-1)^k) * theta.^(2*k) / factorial(2*k);
    end

    % --- pendekatan e^x = sum_{k=0}^{n} x^k / k! ---
    exp_n = 0;
    for k = 0:n
        exp_n = exp_n + x^k / factorial(k);
    end

    h_approx(idx, :) = h0 * (cos_n + exp_n);
end

%% ---------- 4. GALAT ----------
galat_mutlak = abs(repmat(h_eksak, length(n_list), 1) - h_approx);
galat_relatif = galat_mutlak ./ abs(repmat(h_eksak, length(n_list), 1)) * 100;

rerata_mutlak  = mean(galat_mutlak, 2);     % rerata galat mutlak tiap n
rerata_relatif = mean(galat_relatif, 2);    % rerata galat relatif (%) tiap n

%% ---------- 5. TAMPILKAN HASIL ----------
fprintf('==========================================================\n');
fprintf(' n  | Rerata Galat Mutlak | Rerata Galat Relatif (%%) | Layak?\n');
fprintf('==========================================================\n');
for idx = 1:length(n_list)
    if rerata_relatif(idx) < tol
        status = 'LAYAK';
    else
        status = 'TIDAK LAYAK';
    end
    fprintf(' %2d | %19.6f | %24.6f | %s\n', ...
        n_list(idx), rerata_mutlak(idx), rerata_relatif(idx), status);
end
fprintf('==========================================================\n');

%% ---------- 6. GRAFIK 1: PERBANDINGAN 4 FUNGSI ----------
figure(1);
plot(t, h_eksak, 'k-', 'LineWidth', 3); hold on;
plot(t, h_approx(1,:), 'r--', 'LineWidth', 2);
plot(t, h_approx(2,:), 'g-.', 'LineWidth', 2);
plot(t, h_approx(3,:), 'b:',  'LineWidth', 2);
grid on;
xlabel('t'); ylabel('h(t)');
title('Perbandingan Hasil Eksak dan Pendekatan Deret Maclaurin');
legend('f_{eksak}', 'f_{n=23}', 'f_{n=24}', 'f_{n=25}', 'Location', 'northeast');
hold off;

%% ---------- 7. GRAFIK 2: GALAT MUTLAK (agar beda terlihat) ----------
figure(2);
plot(t, galat_mutlak(1,:), 'r-', 'LineWidth', 2); hold on;
plot(t, galat_mutlak(2,:), 'g-', 'LineWidth', 2);
plot(t, galat_mutlak(3,:), 'b-', 'LineWidth', 2);
grid on;
xlabel('t'); ylabel('|h_{eksak} - h_{pendekatan}|');
title('Galat Mutlak terhadap t');
legend('n = 23', 'n = 24', 'n = 25', 'Location', 'east');
hold off;
