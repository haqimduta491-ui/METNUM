% ELIMINASI GAUSS-JORDAN (INTERAKTIF)
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
        Ab(i,:) = Ab(i,:) - m * Ab(k,:);
    end
end
disp('Hasil forward elimination (segitiga atas):');
disp(Ab);

% ---- Backward elimination ----
for k = n:-1:2
    for i = k-1:-1:1
        m = Ab(i,k) / Ab(k,k);
        Ab(i,:) = Ab(i,:) - m * Ab(k,:);
        printf('R%d = R%d - (%g)*R%d\n', i, i, m, k);
    end
end
disp('Setelah backward elimination (matriks diagonal):');
disp(Ab);

% ---- Normalisasi diagonal menjadi 1 ----
for i = 1:n
    Ab(i,:) = Ab(i,:) / Ab(i,i);
end
disp('Bentuk identitas [I|x]:');
disp(Ab);

x = Ab(:,n+1);

disp('=== Solusi ===');
for i = 1:n
    printf('x%d = %.6f\n', i, x(i));
end

disp('Verifikasi (A*x - b, harus mendekati nol):');
disp(A*x - b);
