function [Dips,Dipa,Dops,Dopa,D] = normdet(GEO,MAT,K,W)
%ADIMDET Determinants of the ribbon model as a function of the adimensionnalized 
% wavenumber K = .5*b*k
% frequency W = .5*b*w/vt
[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(GEO,MAT) ;
a = b/h ; % section aspect ratio
x = G./Q ;
alpha = 1 + 1./(x.*xi) ;

% Squared K and W
W2 = W.*W ;
K2 = K.*K ;

% Trigonometric functions 
if isa(W,'math.poly') || isa(K,'math.poly') % polynomial approximations
    N = math.poly(W+K).N ;
    sinc_sqrt = @(X) math.poly.sincX(2*N+4).keeporders(2).comp(X,N+2) ;
    cos_sqrt = @(X) math.poly.cosX(2*N+4).keeporders(2).comp(X,N+2) ;
else
    sinc_sqrt = @(x) sin(sqrt(x))./sqrt(x) ;
    cos_sqrt = @(x) cos(sqrt(x)) ;
end

% IN-PLANE MOTION
% Wavenumbers
Kl2 = x.*W2 - K2 ;
Kt2 = W2 - K2 ;
% Determinants
Dips = 4.*K2.*Kl2.*cos_sqrt(Kt2).*sinc_sqrt(Kl2) + (Kt2-K2).^2.*cos_sqrt(Kl2).*sinc_sqrt(Kt2) ;
Dipa = 4.*K2.*Kt2.*cos_sqrt(Kl2).*sinc_sqrt(Kt2) + (Kt2-K2).^2.*cos_sqrt(Kt2).*sinc_sqrt(Kl2) ;

% OUT-OF-PLANE MOTION
% Wavenumbers
Kst2 = W2 - 3.*xi.*a.^2 - K2 ;
sqrtAlpha2m4Beta2 = sqrt( 12.*xi.*(alpha-1).*a.^2 + W2.*(alpha-2).^2 ) ;
Kb2 = .5.*x.*W.*(alpha.*W + sqrtAlpha2m4Beta2) - K2 ;
Ksa2 = .5.*x.*W.*(alpha.*W - sqrtAlpha2m4Beta2) - K2 ;
% Trigonometric functions
cst = cos_sqrt(Kst2) ; sst = sinc_sqrt(Kst2) ; 
cb = cos_sqrt(Kb2) ; sb = sinc_sqrt(Kb2) ; 
csa = cos_sqrt(Ksa2) ; ssa = sinc_sqrt(Ksa2) ; 
% Determinants
Dops = 4.*x.*K2.*Kb2.*Ksa2.*(Kb2-Ksa2).*sb.*ssa.*cst + sst.*(...
    + Ksa2.*(Kb2+(1-2.*x).*K2).*(2.*x.*K2.*(Kt2-Kst2) + (Kl2-Ksa2).*(Kst2-K2)).*ssa.*cb ...
    - Kb2.*(Ksa2+(1-2.*x).*K2).*(2.*x.*K2.*(Kt2-Kst2) + (Kl2-Kb2).*(Kst2-K2)).*sb.*csa ...
    ) ; % symetric modes
Dopa = 4.*x.*K2.*Kst2.*(Kb2-Ksa2).*sst.*cb.*csa + cst.*( ...
     + (Kb2+(1-2.*x).*K2).*(2.*x.*K2.*(Kt2-Kst2) + (Kl2-Ksa2).*(Kst2-K2)).*sb.*csa ...
     - (Ksa2+(1-2.*x).*K2).*(2.*x.*K2.*(Kt2-Kst2) + (Kl2-Kb2).*(Kst2-K2)).*ssa.*cb ...
     ) ; % antisymetric modes
% Normalize
Dops = Dops./(Kst2+K2) ;
Dopa = Dopa./(Kst2+K2) ;

% Full determinant if needed
if nargout>4 ; D = Dips.*Dipa.*Dops.*Dopa ; end

end

