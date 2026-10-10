% =========================================================
% SOAL 3 : Running time kelima metode dengan tic-toc
%          (hanya algoritma, tanpa input/output)
% Jalankan:  >> soal3
% =========================================================
1;
function x = gauss_elim(A, b)
  % Eliminasi Gauss + pivoting parsial + substitusi mundur
  n = length(b); Ab = [A b];
  for k = 1:n-1
    [~, idx] = max(abs(Ab(k:n, k))); idx = idx + k - 1;
    if idx ~= k, Ab([k idx], :) = Ab([idx k], :); end
    for i = k+1:n
      m = Ab(i,k) / Ab(k,k);
      Ab(i,k:n+1) = Ab(i,k:n+1) - m * Ab(k,k:n+1);
    end
  end
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = (Ab(i,n+1) - Ab(i,i+1:n) * x(i+1:n)) / Ab(i,i);
  end
end

function x = gauss_jordan(A, b)
  % Eliminasi Gauss-Jordan
  n = length(b); Ab = [A b];
  for k = 1:n
    [~, idx] = max(abs(Ab(k:n, k))); idx = idx + k - 1;
    if idx ~= k, Ab([k idx], :) = Ab([idx k], :); end
    Ab(k,:) = Ab(k,:) / Ab(k,k);
    for i = [1:k-1, k+1:n]
      Ab(i,:) = Ab(i,:) - Ab(i,k) * Ab(k,:);
    end
  end
  x = Ab(:, n+1);
end

function x = lu_solve(A, b)
  % Dekomposisi LU (Doolittle): A = L*U, Ly = b, Ux = y
  n = length(b); L = eye(n); U = zeros(n);
  for i = 1:n
    for j = i:n
      U(i,j) = A(i,j) - L(i,1:i-1) * U(1:i-1,j);
    end
    for j = i+1:n
      L(j,i) = (A(j,i) - L(j,1:i-1) * U(1:i-1,i)) / U(i,i);
    end
  end
  y = zeros(n,1);
  for i = 1:n
    y(i) = b(i) - L(i,1:i-1) * y(1:i-1);
  end
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
  end
end

function [x, it, hist] = jacobi(A, b, x0, tol, maxit)
  n = length(b); x = x0; hist = zeros(maxit,n+1);
  for it = 1:maxit
    xn = x;
    for i = 1:n
      idx = [1:i-1, i+1:n];
      xn(i) = (b(i) - A(i,idx) * x(idx)) / A(i,i);
    end
    err = max(abs(xn - x));
    x = xn;
    hist(it,:) = [x' err];
    if err < tol, break; end
  end
  hist = hist(1:it,:);
end

function [x, it, hist] = gauss_seidel(A, b, x0, tol, maxit)
  n = length(b); x = x0; hist = zeros(maxit,n+1);
  for it = 1:maxit
    xo = x;
    for i = 1:n
      idx = [1:i-1, i+1:n];
      x(i) = (b(i) - A(i,idx) * x(idx)) / A(i,i);
    end
    err = max(abs(x - xo));
    hist(it,:) = [x' err];
    if err < tol, break; end
  end
  hist = hist(1:it,:);
end

% ---------------- DATA SPL ----------------
% 5p -  q -  r            = -1
% -p + 5q -  r -  s       =  2
% -p -  q + 4r -  s -  t  =  6
%            r + 4s - 2t  =  2
%       q -  r +  s + 4t  = -1
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b  = [-1; 2; 6; 2; -1];
x0 = zeros(5,1);     % tebakan awal
tol   = 0.001;       % batas toleransi error
maxit = 20;          % maksimum iterasi
nama  = {'p','q','r','s','t'};

N = 1000;   % SPL kecil -> diulang N kali lalu dirata-rata

tic; for k = 1:N, gauss_elim(A, b);                   end; t_gauss = toc / N;
tic; for k = 1:N, gauss_jordan(A, b);                 end; t_gj    = toc / N;
tic; for k = 1:N, lu_solve(A, b);                     end; t_lu    = toc / N;
tic; for k = 1:N, jacobi(A, b, x0, tol, maxit);       end; t_jac   = toc / N;
tic; for k = 1:N, gauss_seidel(A, b, x0, tol, maxit); end; t_gs    = toc / N;

disp('================ SOAL 3 : RUNNING TIME (rata-rata 1000x) ================');
fprintf('%-16s : %.6e detik\n', 'Gauss',          t_gauss);
fprintf('%-16s : %.6e detik\n', 'Gauss-Jordan',   t_gj);
fprintf('%-16s : %.6e detik\n', 'Dekomposisi LU', t_lu);
fprintf('%-16s : %.6e detik\n', 'Jacobi',         t_jac);
fprintf('%-16s : %.6e detik\n', 'Gauss-Seidel',   t_gs);

waktu  = [t_gauss t_gj t_lu t_jac t_gs];
metode = {'Gauss','Gauss-Jordan','LU','Jacobi','Gauss-Seidel'};
[~, urut] = sort(waktu);
fprintf('\nUrutan tercepat -> terlambat: %s\n', strjoin(metode(urut), ' < '));
