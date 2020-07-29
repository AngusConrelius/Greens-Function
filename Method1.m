
k = 20;
m = 5000;
d = 0.1;
B = 1.414/d;
M = m;
nu = 1;
K = 2;
P = 3.14;
p = 2*P/d;
x = 0;
y = 0.05;
n = 0;
z = 0;
X = x - z;
Y = 0.01/d;

G=0;

r = (X^2 + (Y-0*d)^2)^(1/2)
Z=k*r;
H = besselh(nu,K,Z);
G = G + H*exp(1i*0*B*d)
%plot (o,abs(G),'+')
%hold on
%pause (0.01)

for m = 1:M
    r = (X^2 + (Y-m*d)^2)^(1/2); %negitive m = 1 case creats r = 0 there for when used in bessel function creates a nan number
    Z=k*r;
    
    %if z ~= 0
       %H = besselh(nu,K,Z); %avoids the z=0 casec that created NaN
    %end
    
    G = G + H*exp(1i*m*B*d)
    %plot (m,abs(G),'+')
    %hold on
    %pause (0.01)
    r = (X^2 + (Y+m*d)^2)^(1/2);
    Z=k*r;
    H = besselh(nu,K,Z);
    G = G + H*exp(1i*m*B*d)
    %plot (m,abs(G),'+')
    %hold on
    %pause (0.01)
end


Q = vpa(-1i/4)*G %vpa increases precision to 32 decimals instead of base 16

