%% ============================================================
%  RESPONSI 1 PRAKTIKUM METODE NUMERIK
%  Kasus 2 : SPL (Harga Bahan Makanan per Kilogram)
%% ============================================================
clc; clear; close all;

%% ---------- 1. KONVERSI SATUAN KE KG ----------
% Kolom: Beras, Telur, Tepung, Ayam | Baris: Adel, Briza, Caca, Dinda
% 1 hg = 0.1 kg ; 1 dag = 0.01 kg ; 1 g = 1/1000 kg
A = [2.5,       1500/1000, 5,       500/1000;
     40*0.1,    0.75,      1000/1000, 2;
     2000/1000, 15*0.1,    50*0.01, 6000/1000;
     1,         4250/1000, 20*0.1,  9*0.1];
b = [161000; 149000; 251500; 173000];
nama = {'Beras', 'Telur', 'Tepung', 'Ayam'};

%% ---------- 2. SPL DAN MATRIKS ----------
disp('Matriks koefisien A (kg):'); disp(A);
disp('Vektor konstanta b (Rp):'); disp(b);

%% ---------- 3. CEK DOMINAN DIAGONAL (URUTAN ASLI) ----------
n = length(b);
fprintf('\nCek dominan diagonal (urutan asli):\n');
for i = 1:n
    diag_i = abs(A(i,i));
    sisa   = sum(abs(A(i,:))) - diag_i;
    fprintf('Baris %d : |a_ii| = %.2f , jumlah lainnya = %.2f\n', i, diag_i, sisa);
end

%% ---------- 4. PENUKARAN BARIS ----------
urut = [2 4 1 3];            % Briza, Dinda, Adel, Caca
A2 = A(urut, :);
b2 = b(urut);
fprintf('\nCek dominan diagonal (setelah penukaran baris):\n');
dominan = true;
for i = 1:n
    diag_i = abs(A2(i,i));
    sisa   = sum(abs(A2(i,:))) - diag_i;
    fprintf('Baris %d : |a_ii| = %.2f , jumlah lainnya = %.2f\n', i, diag_i, sisa);
    if diag_i <= sisa
        dominan = false;
    end
end
if dominan
    disp('Dominan diagonal terpenuhi -> Gauss-Seidel konvergen');
else
    disp('Dominan diagonal tidak terpenuhi');
end

%% ---------- 5. METODE GAUSS-SEIDEL ----------
x    = zeros(n,1);           % tebakan awal
tol  = 1e-6;
maks = 500;
fprintf('\nIterasi Gauss-Seidel:\n');
fprintf('Iter |    Beras      Telur     Tepung       Ayam | Galat maks\n');
for iter = 1:maks
    x_lama = x;
    for i = 1:n
        jumlah = A2(i,:)*x - A2(i,i)*x(i);
        x(i) = (b2(i) - jumlah) / A2(i,i);
    end
    galat = max(abs(x - x_lama));
    fprintf('%4d | %10.4f %10.4f %10.4f %10.4f | %.6e\n', iter, x, galat);
    if galat < tol
        break;
    end
end

%% ---------- 6. VERIFIKASI DENGAN ELIMINASI GAUSS ----------
M = [A b];
for k = 1:n-1
    [~, p] = max(abs(M(k:n,k)));  p = p + k - 1;
    M([k p], :) = M([p k], :);
    for i = k+1:n
        f = M(i,k) / M(k,k);
        M(i,:) = M(i,:) - f*M(k,:);
    end
end
xg = zeros(n,1);
for i = n:-1:1
    xg(i) = (M(i,n+1) - M(i,i+1:n)*xg(i+1:n)) / M(i,i);
end

%% ---------- 7. HASIL AKHIR ----------
fprintf('\nHasil harga per kilogram:\n');
fprintf('Bahan  | Gauss-Seidel   | Gauss (verifikasi)\n');
for i = 1:n
    fprintf('%-6s | Rp %10.2f | Rp %10.2f\n', nama{i}, x(i), xg(i));
end
fprintf('Jumlah iterasi Gauss-Seidel : %d\n', iter);
fprintf('Cek residu norm(A*x - b)    : %.6e\n', norm(A*x - b));
