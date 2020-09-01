
function [S0]= ZeroSum(l)

m = 7;
k = 20;
d = 0.1;
B = sqrt(2)/d;
M = m;
c = 0.5772157; % PDF page 11
C = 1.2020569; % PDF page 11
p = (2*pi)/d;   % PDF page 4
Y0 = (B^2 - k^2)^(1/2); % Zero Ym
Ym = 0;

Sum2 = 0; % the sum number is related to which sum this is in the entire lattice sum method so they can be combined into one file in the future
for m = -M:M
    if (m == 0)
        
        Sum2 = Sum2;
    else
        Bm = B + m*p; % PDF page 5
        Ym = (Bm^2 - k^2)^(1/2); % PDF page 6
        Sum2 = Sum2 + 1/Ym - 1/p*abs(m) - (k^2 + 2*(B^2))/(2*(p^3)*abs(m)^2);
    end
end



S0 = (-1)-(2i/pi)*(c+log(k/2*p))-(2i/Y0*d)-((2i)*(k^2 + 2*(B^2))/p^3 * d)*C - (2i/d) * Sum2;

end
