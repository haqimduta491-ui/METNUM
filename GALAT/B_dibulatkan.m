clc;
clear;

n = 20;
jumlah_bulat = 0;

for i = 1:n
    hasil_bagi = round((1/i) * 10000) / 10000;
    jumlah_bulat = jumlah_bulat + hasil_bagi;
end

fprintf('Hasil dibulatkan = %.10f\n', jumlah_bulat);
