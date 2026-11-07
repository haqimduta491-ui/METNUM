
%        Python: return (x**2 - 1) / y
function r = f(x, y)
  r = (4*x^3 - 2*x + 1) / (2 + y);
end


function r = solusi_eksak(x)
  r = -2 + sqrt(2*x^4 - 2*x^2 + 2*x + 9);
end


function [tabel, y] = metode_euler(x_awal, y_awal, x_target, N)
  h = (x_target - x_awal) / N;
  x = x_awal;  y = y_awal;
  tabel = zeros(N, 5);                     % [TAMBAH]
  for i = 0:N-1                            % [UBAH] range(N) -> 0:N-1
    fi = f(x, y);
    y_next = y + h * fi;
    tabel(i+1, :) = [i, x, y, fi, y_next]; % [UBAH] tabel.append(...) -> isi baris matriks
    x = x + h;
    y = y_next;
  end
end

% [UBAH] print(f"...") -> printf('...')
function cetak_tabel(t)
  printf('%-3s%-10s%-12s%-14s%-12s\n', 'i', 'xi', 'yi', 'f(xi,yi)', 'yi+1');
  printf('%s\n', repmat('.', 1, 51));
  for r = 1:size(t, 1)
    printf('%-3d%-10.4f%-12.4f%-14.4f%-12.4f\n', t(r, :));
  end
end


x_awal = 0;  y_awal = 1;
x_target = 2;
N = 15;

[tabel, y_euler] = metode_euler(x_awal, y_awal, x_target, N);

printf('     Tabel Iterasi Metode Euler  \n\n');
cetak_tabel(tabel);

y_eksak = solusi_eksak(x_target);
galat_absolut = abs(y_eksak - y_euler);
galat_relatif = galat_absolut / y_eksak * 100;

% [UBAH] teks "sqrt(73)" dihapus karena itu jawaban soal lama; sekarang -2+sqrt(37)
printf('\n       Perbandingan Solusi Eksak\n');
printf('Solusi eksak y(%d) = -2 + sqrt(37) = %.6f\n', x_target, y_eksak);
printf('Metode Euler y(%d)                 = %.6f\n', x_target, y_euler);
printf('Galat absolut                      = %.6f\n', galat_absolut);
printf('Galat relatif                      = %.4f%%\n', galat_relatif);

