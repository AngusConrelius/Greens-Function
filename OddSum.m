function [s3L]= OddSum(l)

d = 2/20;
B = sqrt(2)/d;
m = 7;
M = m;
p = 2*pi/d;
k = 20;
z = zeta(2*l + 1);
Bm = B + m*p;
Ym = 0;
nYm = 0;
Y0 = (B^2 - k^2)^(1/2);
theta0 = sin(1)^(-1) * (0/k);


Term12 = 2i*(exp(-1i*(2*l-1)*theta0)) / Y0*d

Sum5 = 0  % the sum number is related to which sum this is in the entire lattice sum method so they can be combined into one file in the future
for m = 1:M
    Bm = B + m*p;
    theta = sin(1)^(-1) * (Bm/k);
    Ntheta = sin(1)^(-1) * (-Bm/k);
    nYm = (-Bm^2 - k^2)^(1/2);
    Ym = (Bm^2 - k^2)^(1/2);
    
    Sum5 = Sum5 + (exp(-1i*(2*l-1)*theta) / (Ym*d)) - (exp(-1i*(2*l-1)*Ntheta) / (nYm*d)) + ((1i*((-1)^l)*B*d*l) / ((m*pi)^2) * (k/2*m*p)^(2*l-1));
end

Term13 = (2*(-1)^(l) * B*d*l / pi^2)*(k/2*p)^2*l-1 * z

Sum6 = 0;
for m = 0:(l-1)
    
    Term14 = (-1)^(m) * 2^(2*m)
    
    format long
    fact4 = factorial(l+m-1)
    
    
    format long
    fact5 = factorial(2*m+1)
    
    
    format long
    fact6 = factorial(l-m-1)
    
    Ber = bernoulli(2*m + 1, B/p)
    Sum6 = Sum6 + (Term14*fact4)/(fact5*fact6) * (p/k)^2*m+1 * Ber;
    
end

s3L = Term12 + (2i*Sum5) + Term13 - (2/pi) * Sum6

end






