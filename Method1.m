
k = 20;
m = 5000;
d = 2/20;

B = sqrt(2)/d;   

M = m;
nu = 0;         
K = 2;

p = 2*pi/d;

X = 0;
Y = 0.01*d;
G = 0;


r = sqrt((X^2 + (Y-0*d)^2));
Z = k*r;
G = G + (-1i/4)*besselh(nu,K,Z)*exp(1i*0*B*d);
%plot (o,abs(G),'+')
%hold on
%pause (0.01)

S(1) = G;

for m = 1:M
    
    
    r = (X^2 + (Y-m*d)^2)^(1/2);
    Z=k*r;
    G = G + (-1i/4)*besselh(nu,K,Z)*exp(1i*m*B*d);
    %plot (m,abs(G),'+')
    %hold on
    %pause (0.01)
    
    
    r = (X^2 + (Y-(-m)*d)^2)^(1/2);
    Z=k*r;
   
    G = G + (-1i/4)*besselh(nu,K,Z)*exp(1i*(-m)*B*d);
    %plot (m,abs(G),'+')
    %hold on
    %pause (0.01)
    
    S(m+1) = G;
    
   % if  mod(m,100) == 0
        
      %  plot(abs(S)); hold on
        
      %  plot(real(S));
        
       % plot(imag(S));
        
      %  hold off
        
        
        
      %  pause (0.01)
        
   % end
end

G



