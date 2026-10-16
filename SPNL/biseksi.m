clc; clear all; close all;

% ======================== INPUT (UBAH DI SINI) =======================
pilih = 3;        % pilih persamaan: 1, 2, atau 3
n     = 20;       % maksimal iterasi
tol   = 0.0001;   % batas error toleransi

% Daftar persamaan beserta interval awal [a, b]
% (untuk persamaan lain, ubah salah satu case di bawah)
switch pilih
  case 1   % y = 3x^3 - 5x^2 + 4x - 7 dengan y = 1
    f = @(x) 3*x.^3 - 5*x.^2 + 4*x - 7 - 1;
    a = 1;   b = 2;
  case 2   % y = e^x - x^2 dengan y = 0
    f = @(x) exp(x) - x.^2;
    a = -1;  b = 0;
  case 3   % y = log(x+2) + 2x - 3 dengan y = 0 (log basis 10)
    f = @(x) log10(x+2) + 2*x - 3;
    a = 1;   b = 2;
end
% =====================================================================

% ------------------------- PROSES BISEKSI ----------------------------
fa = f(a);
fb = f(b);

if fa*fb > 0
  error('Warning: f(a) dan f(b) sama tanda. Ganti interval a dan b.');
end

fprintf('\n=== METODE BISEKSI | Persamaan %d | interval [%g, %g] ===\n', pilih, a, b);
fprintf('Iter\ta\t\tb\t\tfa\t\tfb\t\tc\t\tfc\n');

for i = 1:n
  c  = (a+b)/2;          % titik tengah
  fc = f(c);
  fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\n', i, a, b, fa, fb, c, fc);

  % berhenti jika error mencapai toleransi
  if abs(c-a) <= tol || abs(c-b) <= tol
    break;
  end

  if fa*fc == 0          % akar tepat ditemukan
    break;
  elseif fa*fc < 0       % akar di antara a dan c
    b = c;  fb = f(b);
  else                   % akar di antara c dan b
    a = c;  fa = f(a);
  end
end

fprintf('\nAkar hampiran = %.6f\n', c);
fprintf('Jumlah iterasi = %d\n', i);
