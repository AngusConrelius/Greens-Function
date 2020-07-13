
k = 2;
m=5000;
B = 1.414;
d=1;
M = m; 
nu = 2;
K=2;
Z=2;
x = 1;
y = 1;
n = 1;
X = x - Z;
Y = y - n; % n is a metric tensor

z = zeta(sym([0.7 i 4 11/3]));

r = (X^2 + (Y-m*d^2))^(1/2);

H = besselh(nu,K,Z);

for m = -M:M
   G = H*(k*r)*exp(1)^(1i*m*B*d);
end

G = (-1/4)*G

