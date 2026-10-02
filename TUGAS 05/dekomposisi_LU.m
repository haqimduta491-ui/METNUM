% DEKOMPOSISI LU (L diagonal = 1), INTERAKTIF
clc; clear;

n = input('Masukkan jumlah persamaan (n): ');

A = zeros(n,n);
b = zeros(n,1);

disp('Masukkan matriks koefisien A:');
for i = 1:n
    for j = 1:n
        A(i,j) = input(sprintf('A(%d,%d) = ', i, j));
    end
end

disp('Masukkan vektor b:');
for i = 1:n
    b(i) = input(sprintf('b(%d) = ', i));
end

U = A;
L = eye(n);

% ---- Forward elimination, simpan pengali m ke L ----
for k = 1:n-1
    if U(k,k) == 0
        error('Pivot nol pada baris %d, program berhenti (butuh pivoting).', k);
    end
    for i = k+1:n
        m = U(i,k) / U(k,k);
        L(i,k) = m;
        U(i,:) = U(i,:) - m * U(k,:);
        printf('m%d%d = %g\n', i, k, m);
    end
end

disp('Matriks L:');
disp(L);
disp('Matriks U:');
disp(U);

disp('Cek L*U:');
disp(L*U);
disp('A:');
disp(A);
if norm(A - L*U) < 1e-9
    disp('Terbukti: A = L*U');
else
    disp('Peringatan: A != L*U (cek input / pivot nol)');
end

% ---- Substitusi maju: L y = b ----
y = zeros(n,1);
for i = 1:n
    y(i) = (b(i) - L(i,1:i-1) * y(1:i-1)) / L(i,i);
end
disp('Hasil substitusi maju, y:');
disp(y);

% ---- Substitusi mundur: U x = y ----
x = zeros(n,1);
for i = n:-1:1
    x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
end

disp('=== Solusi ===');
for i = 1:n
    printf('x%d = %.6f\n', i, x(i));
end

disp('Verifikasi (A*x - b, harus mendekati nol):');
disp(A*x - b);
