function [kips,kipa,kops,kopa] = cutoffK(n,GEO,MAT)
% Compute the n first cutoff wavenumbers k related to a zero frequency w=0
[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(GEO,MAT) ;

% IN-PLANE MOTION: sh(k.b)+m.k.b = 0 
% Find the roots with the polynomial approximation
nn = 0:4*n ; % degree of aproximation
% Approximate sh(x) with the series x^(2n+1)/factorial(2n+1)
p = sparse(2*nn+2,1,1./factorial(2*nn+1)) ;
% Symmetric motion (m=1)
    ps = p + sparse(2,1,1,numel(p),1) ;
    kips = roots(flip(ps))/b ;
    kips = sort(kips,'ascend','comparisonmethod','abs') ;
    kips = kips(1:n) ;
% Anti-symmetric motion (m=-1)
    pa = p - sparse(2,1,1,numel(p),1) ;
    kipa = roots(flip(pa))/b ;
    kipa = sort(kipa,'ascend','comparisonmethod','abs') ;
    kipa = kipa(1:n) ;
    
return
%% OUT-OF-PLANE MOTION
    % quantities of interest
    xi = pi^2/12 ;
    vinf = sqrt(xi).*vt ;
    w0 = sqrt(12)/GEO.h*vinf ;
    % wavenumbers (isotropic part mu_i² = k_i² + k²)
    must2 = (w.^2-w0.^2)./vt.^2 ; % transverse shear
    alpha = (vl./vinf).^2 + 1 ;
    beta = 4.*(1-(w0./w).^2).*(vl./vinf).^2 ;
    gamma = sqrt(alpha.^2 - beta) ;
    mub2 = .5*(w./vl).^2.*(alpha + gamma) ; % bending
    musa2 = .5*(w./vl).^2.*(alpha - gamma) ; % axial shear
    % width-direction wavenumbers & trigo
    kb = sqrt(mub2) ; cb = cos(GEO.b./2.*kb) ; sb = sin(GEO.b./2.*kb) ; 
    ksa = sqrt(musa2) ; csa = cos(GEO.b./2.*ksa) ; ssa = sin(GEO.b./2.*ksa) ;
    kst = sqrt(must2) ; cst = cos(GEO.b./2.*kst) ; sst = sin(GEO.b./2.*kst) ;
    % precomputation
    Db = ... Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
         ... kb.*(alpha - 2 - gamma) ...
         sqrt(alpha + gamma).*(alpha - 2 - gamma) ...
         ; 
    Dsa = ... Ssa3.*(Ab1.*Sst2-Ab2.*Sst1) ...
          ... ksa.*(alpha - 2 + gamma) ...
          sqrt(alpha - gamma).*(alpha - 2 + gamma) ...
          ; 
    DbpDsa = ... (kb+ksa).*(alpha-2) - sqrt(alpha.^2 - beta).*(kb-ksa) ...
             sqrt(alpha+sqrt(beta)).*(alpha-2) - sqrt((alpha.^2 - beta).*(alpha-sqrt(beta))) ...
             ; % Db+Dsa
    DbmDsa = ... (kb-ksa).*(alpha-2) - sqrt(alpha.^2 - beta).*(kb+ksa) ...
             sqrt(alpha-sqrt(beta)).*(alpha-2) - sqrt(alpha.^2 - beta).*sqrt(alpha+sqrt(beta)) ...
             ; % Db-Dsa
    % determinants
%     detBs = sst.*(cb.*ssa.*Db - csa.*sb.*Dsa) ; % symmetric modes
%     detBa = cst.*(sb.*csa.*Db - ssa.*cb.*Dsa) ; % antisymmetric modes
    detBs = sst.*(DbmDsa.*sin(GEO.b./2.*(kb+ksa)) - DbpDsa.*sin(GEO.b./2.*(kb-ksa))) ; % symmetric modes
    detBa = cst.*(DbmDsa.*sin(GEO.b./2.*(kb+ksa)) + DbpDsa.*sin(GEO.b./2.*(kb-ksa))) ; % antisymmetric modes
    % Analytic cases
    n = 0:7 ;
    wst0 = sqrt(w0.^2 + (vt.*n.*pi/GEO.b).^2) ; % {sin|cos}(b/2*must)=0

clf ; 
plot(w,-log(abs(detAs))) ; 
plot(w,-log(abs(detAa))) ;
plot(wl0.*[1;1],get(gca,'ylim')','-')
plot(wt0.*[1;1],get(gca,'ylim')','-')

clf ; 
plot(w,-log(abs(detBs))) ; 
plot(w,-log(abs(detBa))) ;
plot(w,-log(abs((kb-ksa) - 7*pi/GEO.b)))
plot(wst0.*[1;1],get(gca,'ylim')',':')
% plot(wt0.*[1;1],get(gca,'ylim')','-')
    
end

