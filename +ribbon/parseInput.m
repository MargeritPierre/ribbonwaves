function [h,b,Q,G,nu,rho,xi] = parseInput(GEO,MAT)

h = GEO.h ;
b = GEO.b ;

Q = MAT.Q ;
G = MAT.G ;
nu = 1-2.*G./Q ;
rho = MAT.rho ;

if isfield(MAT,'xi') 
    xi = MAT.xi ;
elseif 0 % matching the first cutoff frequency
% lmbda0 = 2*h = ct/f0 =  2*pi*ct/w0
% -> w0 = pi*ct/h = sqrt(xi)*ct*sqrt(12)/h
    xi = pi^2/12 ;
else % matching the asymptotic wave velocity (rayleigh ?)
% vR = vinf = sqrt(xi)*ct
    xi = analytical.rayleigh(nu).^2 ;
end

end