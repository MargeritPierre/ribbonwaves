function [kips,kipa,kops,kopa] = computeK(GEO,MAT,w)
% COMPUTEW Computes the wavenumbers k associated to a ribbon and a series of
% frequencies w
[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(GEO,MAT) ;
a = b/h ; % section aspect ratio
x = G./Q ;
alpha = 1 + 1./(x.*xi) ;

% Estimate the polynomial expansion order to reach convergence
Nrange = [3 80] ; % range of approximation orders
tol = sqrt(eps) ; % desired tolerance
W2max = 3*xi*a.^2 ; % maximum normalized frequency as height-cutoff frequency
Kb2max = .5*x.*(W2max.*alpha + sqrt((alpha-2).^2.*W2max.^2 + 12*xi.*a.^2.*W2max.*(alpha-1))) ; % maximum bending wavenumber
K2max = abs(Kb2max) ;
fun = @(N)log(gamma(2*N+2)) - N*log(K2max) + log(tol) ; % convergence of the sinc function
N = 2*ceil(fzero(fun,Nrange)) ;

% Compute the polynomial expansion of the determinants at order N
vt = sqrt(G./rho) ;
W = .5.*b.*w./vt ;
K = math.poly.X(N) ;
[Dips,Dipa,Dops,Dopa] = ribbon.normdet(GEO,MAT,K,W) ;

% Keep only specific orders of K^2 (others vanish)
Dips.coeffs = Dips.coeffs(1:2:2*N+1,:) ;
Dipa.coeffs = Dipa.coeffs(1:2:2*N+1,:) ;
Dops.coeffs = Dops.coeffs(1:2:2*N+1,:) ;
Dopa.coeffs = Dopa.coeffs(1:2:2*N+1,:) ;
% Compute the roots of the determinants
Kips2 = roots(Dips) ;
Kipa2 = roots(Dipa) ;
Kops2 = roots(Dops) ;
Kopa2 = roots(Dopa) ;
% Keep only converged values
Kips2(abs(Kips2)>K2max) = NaN ; 
Kipa2(abs(Kipa2)>K2max) = NaN ; 
Kops2(abs(Kops2)>K2max) = NaN ; 
Kopa2(abs(Kopa2)>K2max) = NaN ; 
% Convert to frequency
kips = (2/b)*sqrt(Kips2) ;
kipa = (2/b)*sqrt(Kipa2) ;
kops = (2/b)*sqrt(Kops2) ;
kopa = (2/b)*sqrt(Kopa2) ;

end

