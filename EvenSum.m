
function [s2L]= EvenSum(l)

d = 3;
M = 7;
m = M;
p = (2*pi)/d;
B = sqrt(2)/d;
k = 20;

d = 2/20;

Bm = 0;
Ym = 0;
nYm = 0;
Y0 = (B^2 - k^2)^(1/2);
theta0 = sin(1)^(-1) * (0/k);


Sum3 = 0; % the sum number is related to which sum this is in the entire lattice sum method so they can be combined into one file in the future
for m = 1:M
    theta = sin(1)^(-1) * (Bm/k); % PDF page 11
    Ntheta = sin(1)^(-1) * (-Bm/k); % Negative theta
    Bm = B + m*p; %PDF page 5
    nYm = (-Bm^2 - k^2)^(1/2); % Negative Ym
    Ym = (Bm^2 - k^2)^(1/2);
    
    Sum3 = Sum3 + (exp(-2i*l*theta) / Ym*d) + (exp(2i*l*Ntheta) / nYm*d) - (((-1)^l / m*pi) * ((k / 2*m*p)^2*l))
end

z  = zeta(2*l + 1); % PDF page 11

Term3 = (2i*exp(-2i*l*theta0) / Y0*d)

Term4 = ((2i*(-1)^l) / pi) * ((k/2*p)^2*l) * z + (1i/l*pi)



Sum4 = 0;
for m = 1:l
    
    Term5 = (l+m-1)
    format long
    fact = factorial(Term5)
    
    Term6 = ((-1)^m)*(2^(2*m))*fact
    
    Term7 = (2*m);
    format long
    fact2 = factorial(Term7)
    
    Term8 = (l-m);
    format long
    fact3 = factorial(Term8)
    
    Term9 = fact2 * fact3
    Term10 = Term6/Term9
    
    
    Ber = bernoulli(2*m, B/p) % PDF page 16
    
    Term11 = ((p/k)^(2*m))*Ber;
    Sum4 = Sum4 + (Term10 * Term11)
    
end




s2L = Term3 - 2i*Sum3 - Term4 + (1i/pi)*Sum4

end



