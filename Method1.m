
k = 20;
m=5000;
d=0.1;
B = 1.414/d;
M = m;
nu = 1;
K=2;
P=3.14
p=2*P/d
x = 0;
y = 0.1;
n = 0;
z = 0;
X = x - z;
Y = y - n; 

G=0;

r = (X^2 + (Y-0*d^2))^(1/2);
Z=k*r;
H = besselh(nu,K,Z);
G = G + H*exp(1)^(1i*m*B*d);

for m = 1:M
    r = (X^2 + (Y-m*d^2))^(1/2);
    Z=k*r;
    H = besselh(nu,K,Z);
    G = G + H*exp(1)^(1i*m*B*d);
    %   plot (m,abs(G),'+')p
    %   hold on
    %  pause (0.01)
    r = (X^2 + (Y+m*d^2))^(1/2);
    Z=k*r;
    H = besselh(nu,K,Z);
    G = G + H*exp(1)^(1i*m*B*d);
    % plot (m,abs(G),'+')
    %  hold on
    %  pause (0.01)
end


Q = (-1i/4)*G

