clear; clc;

f = @(x, y) (4*x^3 - 2*x + 1) / (2 + y);

x_awal = 0;  y_awal = 1;
x_target = 2;
N = 15;
h = (x_target - x_awal) / N;

% ----- Iterasi RK4 -----
x = x_awal;  y = y_awal;
tabel = zeros(N, 8);
for i = 0:N-1
  k1 = f(x, y);
  k2 = f(x + h/2, y + k1*h/2);
  k3 = f(x + h/2, y + k2*h/2);
  k4 = f(x + h,   y + k3*h);
  y_next = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
  tabel(i+1, :) = [i, x, y, k1, k2, k3, k4, y_next];
  x = x + h;
  y = y_next;
end
y_rk4 = y;

% ----- Cetak tabel -----
printf('     Tabel Iterasi Metode Runge-Kutta Orde-4\n\n');
printf('%-3s%-9s%-10s%-9s%-9s%-9s%-9s%-10s\n', 'i', 'xi', 'yi', 'k1', 'k2', 'k3', 'k4', 'yi+1');
printf('%s\n', repmat('.', 1, 68));
for r = 1:N
  printf('%-3d%-9.4f%-10.4f%-9.4f%-9.4f%-9.4f%-9.4f%-10.4f\n', tabel(r, :));
end

% ----- Solusi eksak & galat -----
y_eksak = -2 + sqrt(2*x_target^4 - 2*x_target^2 + 2*x_target + 9);
galat_absolut = abs(y_eksak - y_rk4);
galat_relatif = galat_absolut / y_eksak * 100;

printf('\n       Perbandingan Solusi Eksak\n');
printf('Solusi eksak y(2) = -2 + sqrt(37) = %.6f\n', y_eksak);
printf('Metode RK4   y(2)                 = %.6f\n', y_rk4);
printf('Galat absolut                     = %.6f\n', galat_absolut);
printf('Galat relatif                     = %.6f %%\n', galat_relatif);

