clear; clc;

f = @(x, y) (4*x^3 - 2*x + 1) / (2 + y);

x_awal = 0;  y_awal = 1;
x_target = 2;
N = 15;
h = (x_target - x_awal) / N;

% ----- Iterasi Heun -----
x = x_awal;  y = y_awal;
tabel = zeros(N, 6);
for i = 0:N-1
  k1 = f(x, y);
  k2 = f(x + h, y + k1*h);
  y_next = y + (h/2) * (k1 + k2);
  tabel(i+1, :) = [i, x, y, k1, k2, y_next];
  x = x + h;
  y = y_next;
end
y_heun = y;

% ----- Cetak tabel -----
printf('     Tabel Iterasi Metode Heun (RK-2)\n\n');
printf('%-3s%-9s%-11s%-11s%-11s%-10s\n', 'i', 'xi', 'yi', 'k1', 'k2', 'yi+1');
printf('%s\n', repmat('.', 1, 55));
for r = 1:N
  printf('%-3d%-9.4f%-11.4f%-11.4f%-11.4f%-10.4f\n', tabel(r, :));
end

% ----- Solusi eksak & galat -----
y_eksak = -2 + sqrt(2*x_target^4 - 2*x_target^2 + 2*x_target + 9);
galat_absolut = abs(y_eksak - y_heun);
galat_relatif = galat_absolut / y_eksak * 100;

printf('\n       Perbandingan Solusi Eksak\n');
printf('Solusi eksak y(2) = -2 + sqrt(37) = %.6f\n', y_eksak);
printf('Metode Heun  y(2)                 = %.6f\n', y_heun);
printf('Galat absolut                     = %.6f\n', galat_absolut);
printf('Galat relatif                     = %.4f %%\n', galat_relatif);

