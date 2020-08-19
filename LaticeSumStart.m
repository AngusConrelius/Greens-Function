EL = 2;
E = 1;

SL = 1;
l = 4;
O = 1;
o = 1;
m = 100;
y = 0.1;
p = 1;
k = 10;
B = 11;
b = 2;
C = 1.2020569;
i = 1;
c = 1;
d = 0.1;
K = 2;
M = m;
X = 0; %X = r cos θ
Y = 0.1; % Y = r sin θ.
nu = 1;
L = 4;
P = 3.14;
r = X^2 + Y^2;

Sum1 = 0;
for l=0:L
    
    %JL = -1i*(P*l)^-1 * (r/m*d)^l % current place holder untill right JL is found 
    
    for m = 1:5000
        Z=K*m*d;
        H = besselh(nu,K,Z);
        SL = H*(exp(1i*m*B*d) + ((-1)^l) * exp(1)^(-1i*B*m*d))
    end
    
    Sum1 = EL*SL*JL*(k*r)*cos(l)*(P/2 - 0);
    %plot (m,abs(G),'+')
    %hold on
    %pause (0.01)
end
Z=K*r;
H = besselh(nu,K,Z);
Term = -i/4 *(H + Sum1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Term2 = (-1)-(2i/pi)*(c+log(k/2*p))-(2*i/y*d)-((2*i)*(k^2 + B^2)/p*3*d)*C - (2*i/d) %not summed

Sum2 = 0;
for m = -M:M
    Sum2 =(1/y*m)-(1/p*abs(m))-(K^2 + 2*b^(2))/((2*p^3)*abs(m)) % the fuction that is summed
end



S0 = Term2 * Sum2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = zeta(2,3)

Term3 = (2i*exp(-2i*l*O)/y*d)

Sum3 = 0;
for m = 1:M
    Sum3 = (exp(1)^-2*i*l*o*m /y*m*d ) - (exp(1)^i*2*l-1*o-m / y-m*d) - (-1)^l/m*pi*((k/2*m*p)^2*l);
end

Term4 = (((2i*(-1)^l) / P)*(k/2*p)^2*l) * z + i/l*P



Sum4 = 0;
for m = 1:l
    
    Ber = bernoulli(0);
    
    Term5 = (l+m-1);
    format long
    Fact = factorial(Term5)
    
    Term7 = (2*m);
    format long
    Fact2 = factorial(Term7)
    
    Term8 = (l-m);
    format long
    Fact3 = factorial(Term8)
    
    Term6 = ((-1)^m)*(2^(2*m))*Fact
    
    Term9 = Fact2 * Fact3;
    
    Term10 = Term6/Term9;
    
    Term11 = ((p/k)^(2*m))*Ber*(B/p);
    
    Sum4 = Term10 * Term11;
end

S2L = Term3 - 2i*Sum3 - Term4 + 1i*Sum4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = zeta(2*l-1);

Term12 = 2i*(exp(1)^-2i*l*O) / Y*d

Sum5 = 0;
for m = 1:M
    Sum5 = exp(1)^(-1i*2*l-1)*O / (Y*d) - exp(1)^(-1i*2*l-1*o) / (y*d) + ((1i*((-1)^l)*B*d*l) / (m*P)^2) * (k/2*m*p)^(2*l-1)
end

Term13 = ((-1)^(l) *B*d*l / P^2)*(k/2*p)^2*l-1 * z;

Ber = bernoulli(0);

Sum6 = 0;
for m = 0:(l-1)
    
    Term14 = (-1)^(m) * 2^2*m
    
    format long
    fact4 = factorial(1+m-1)
    
    
    format long
    fact5 = factorial(2*m+1)
    
    
    format long
    fact6 = factorial(l-m-1)
    
    Ber = bernoulli(2*m+1)
    Sum6 = (Term14*fact4 / fact5*fact6) * (p/k)^2*m+1 * Ber * (B/p);
    
end

S3L = Term12 + 2i*Sum5 + Term13 - (2/P)*Sum6



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       why dose the table only have one result. how are S0, S2L, and S3L
%       working together
%
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











