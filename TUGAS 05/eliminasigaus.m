% ELIMINASI GAUSS + SUBSTITUSI MUNDUR (INTERAKTIF)
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

Ab = [A b];
disp('Matriks augmented awal [A|b]:');
disp(Ab);

% ---- Forward elimination ----
for k = 1:n-1
    if Ab(k,k) == 0
        error('Pivot nol pada baris %d, program berhenti (butuh pivoting).', k);
    end
    for i = k+1:n
        m = Ab(i,k) / Ab(k,k);
        Ab(i,k:n+1) = Ab(i,k:n+1) - m * Ab(k,k:n+1);
        printf('m%d%d = %g  ->  R%d = R%d - (%g)*R%d\n', i, k, m, i, i, m, k);
    end
    printf('\nMatriks setelah langkah k = %d:\n', k);
    disp(Ab);
end

% ---- Substitusi mundur ----
x = zeros(n,1);
x(n) = Ab(n,n+1) / Ab(n,n);
for i = n-1:-1:1
    x(i) = (Ab(i,n+1) - Ab(i,i+1:n) * x(i+1:n)) / Ab(i,i);
end

disp('Matriks segitiga atas [U|c]:');
disp(Ab);

disp('=== Solusi ===');
for i = 1:n
    printf('x%d = %.6f\n', i, x(i));
end

disp('Verifikasi (A*x - b, harus mendekati nol):');
disp(A*x - b);
