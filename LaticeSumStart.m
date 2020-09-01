EL = 2;
SL = 1;
m = 7;
d = 2/20;
p = 2*pi/d;
B = sqrt(2)/d;
c = 0.5772157; % in PDF on page 11
C = 1.2020569; % in PDF on page 11
K = 2;
k = 20;
M = m;
X = 0; % X = r cos θ. given in PDF on page 10 but was not implemented
Y = 0.1; % Y = r sin θ. given in PDF on page 10 but was not implemented
nu = 0;
l = 3;
L = 3;
r = X^2 + Y^2;


Sum1 = 0
for l = 0:L
    
    for m = 1:5000 % in the document this is said to be inf. 5000 is a place holder I chose
        Z=k*m*d;
        H = besselh(nu,K,Z);
        SL = SL + (H*(exp(1i*m*B*d) + ((-1)^l) * exp(-1i*B*m*d))) % SL is given on page 10 of PDF
    end
    
    Z=k*r;
    JL = besselj(nu,Z); % was not found directly in document but belived to be a bessel j function
    Bm = B + m*p
    theta = sin(1)^(-1) * (Bm/k);% PDF page 11
    Sum1 = Sum1 + EL*SL*JL*cos(l)*(pi/2 - theta)
end

Z=k*r;
H = besselh(nu,K,Z);
Term = -1i/4 * (H + Sum1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S0 = 0;
s2L = 0;
s3L = 0;
for l = 0:L
    
    if (l==0)
        
        Term = Term + ZeroSum(l)
        
    elseif mod(l,2) == 0
        
        Term = Term + EvenSum(l)
        
    else
        
        Term = Term + OddSum(l)
        
    end
    
end

Term
