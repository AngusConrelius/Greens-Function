
k= 2;
i=1;
m=5000;
B = 1.414;
d=1;
M = m;
nu = 2;
K=2;
Z=1;
X=0;
Y=0;

r = X^2 + ((Y-m*d)^2)^1/2;

H = besselh(nu,K,Z);

for m = -M:M
   G = H*(k*r)*exp(1)^i*m*B*d;
end

U = @(X,Y,M) (-1/4)*G;

fplot(U,[-M M])