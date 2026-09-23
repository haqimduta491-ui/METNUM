# Menghitung nilai sin(x)
x = input('Input nilai x : ');
n = input('Input nilai n : ');

# Perhitungan nilai sin(x) dengan cara eksak
eksak = sin(x);

# Perhitungan nilai sin(x) dengan deret Taylor
p = 0;
for i = 0:n
  p = p + ((-1)^i * x^(2*i+1)) / factorial(2*i+1);
end

galat = abs(eksak - p);

fprintf('Eksak\t\t= %5.15f\n', eksak)
fprintf('Pendekatan\t= %5.15f\n', p)
fprintf('Error\t\t= %5.15f\n', galat)
