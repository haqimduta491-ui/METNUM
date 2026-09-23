clc;
clear;

n = 20;
i_vektor = 1:n;
jumlah_sum = sum(1 ./ i_vektor);

fprintf('Hasil tanpa loop (sum) = %.10f\n', jumlah_sum);
