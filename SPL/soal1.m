% =========================================================
% SOAL 1 : Solusi SPL dengan Gauss, Gauss-Jordan, LU,
%          Jacobi, dan Gauss-Seidel (maks iterasi 20, tol 0.001)
% Jalankan:  >> soal1
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

x_gauss = gauss_elim(A, b);
x_gj    = gauss_jordan(A, b);
x_lu    = lu_solve(A, b);
[x_jac, it_jac, h_jac] = jacobi(A, b, x0, tol, maxit);
[x_gs,  it_gs,  h_gs ] = gauss_seidel(A, b, x0, tol, maxit);

disp('================ SOAL 1 : SOLUSI SPL ================');
fprintf('%-6s %-14s %-14s %-14s %-14s %-14s\n', ...
        'Var','Gauss','Gauss-Jordan','LU','Jacobi','Gauss-Seidel');
for i = 1:5
  fprintf('%-6s %-14.8f %-14.8f %-14.8f %-14.8f %-14.8f\n', nama{i}, ...
          x_gauss(i), x_gj(i), x_lu(i), x_jac(i), x_gs(i));
end

fprintf('\nJacobi       : %d iterasi, selisih terakhir = %.6f\n', it_jac, h_jac(end,end));
fprintf('Gauss-Seidel : %d iterasi, selisih terakhir = %.6f\n', it_gs,  h_gs(end,end));
if h_jac(end,end) >= tol, disp('  (Jacobi belum mencapai toleransi dalam 20 iterasi)'); end
if h_gs(end,end)  >= tol, disp('  (Gauss-Seidel belum mencapai toleransi dalam 20 iterasi)'); end

disp(' ');
disp('Tabel iterasi Jacobi  [iter | p q r s t | selisih maks]:');
disp([(1:it_jac)' h_jac]);
disp('Tabel iterasi Gauss-Seidel  [iter | p q r s t | selisih maks]:');
disp([(1:it_gs)' h_gs]);

fprintf('\nHimpunan penyelesaian (Gauss):\n');
fprintf('  {(p,q,r,s,t)} = {(%.6f, %.6f, %.6f, %.6f, %.6f)}\n', x_gauss);
