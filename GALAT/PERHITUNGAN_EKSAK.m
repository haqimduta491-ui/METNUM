clc;
clear;

n = 20;
jumlah_eksak = 0;

for i = 1:n
    jumlah_eksak = jumlah_eksak + 1/i;
end

fprintf('Hasil eksak = %.10f\n', jumlah_eksak);
