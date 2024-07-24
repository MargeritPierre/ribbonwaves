function [wips,wipa,wops,wopa] = computeW(GEO,MAT,k)
% COMPUTEW Computes the frequencies associated to a ribbon and a series of
% wavenumbers k
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
W = math.poly.X(N) ;
K = .5.*b.*k ;
[Dips,Dipa,Dops,Dopa] = ribbon.normdet(GEO,MAT,K,W) ;

% Keep only specific orders of W^2 (others vanish)
Dips.coeffs = Dips.coeffs(3:2:2*N+1,:) ;
Dipa.coeffs = Dipa.coeffs(3:2:2*N+1,:) ;
Dops.coeffs = Dops.coeffs(2:2:2*N,:) ;
Dopa.coeffs = Dopa.coeffs(2:2:2*N,:) ;
% Compute the roots of the determinants
Wips2 = roots(Dips) ;
Wipa2 = roots(Dipa) ;
Wops2 = roots(Dops) ;
Wopa2 = roots(Dopa) ;
% Keep only converged values
Wips2(abs(Wips2)>W2max) = NaN ; 
Wipa2(abs(Wipa2)>W2max) = NaN ; 
Wops2(abs(Wops2)>W2max) = NaN ; 
Wopa2(abs(Wopa2)>W2max) = NaN ; 
% Convert to frequency
vt2 = G./rho ;
wips = (2/b)*sqrt(vt2.*Wips2) ;
wipa = (2/b)*sqrt(vt2.*Wipa2) ;
wops = (2/b)*sqrt(vt2.*Wops2) ;
wopa = (2/b)*sqrt(vt2.*Wopa2) ;

end

