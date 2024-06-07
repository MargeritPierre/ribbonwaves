%% DETERMINANT COMPUTATION
function [D,detAs,detAa,detBs,detBa] = determinant(h,b,Q,G,rho,w,k,xi,normalize)
if nargin<8 || isempty(xi) ; xi = pi^2/12 ; end % shear correction factor
if nargin<9 || isempty(normalize) ; normalize = false ; end
nu = 1-2.*G./Q ;

% IN-PLANE MOTION
% phase velocities
vl = sqrt(Q./rho) ;
vt = sqrt(G./rho) ;
% wavenumbers
kl = sqrt((w./vl).^2 - k.^2) ; % longitudinal
kt = sqrt((w./vt).^2 - k.^2) ; % transverse
% precompute
cl = cos(b./2.*kl) ; sl = sin(b./2.*kl) ; 
ct = cos(b./2.*kt) ; st = sin(b./2.*kt) ; 
Al = 4*(k.^2.*kl.*kt) ; At = (kt.^2-k.^2).^2 ;
% determinants
detAs = Al.*sl.*ct + At.*st.*cl ; % symmetrical modes
detAa = Al.*st.*cl + At.*sl.*ct ; % antisymmetrical modes
% normalized
if normalize
    detAs = detAs./Al...
            .*abs(k).^2.*abs(w).^-2 ...
            .*exp(b/2*imag(kt+kl)) ...
            ; % symmetrical modes
    detAa = detAa./Al...
            .*abs(k).^2.*abs(w).^-2 ...
            .*abs(kt).^1.*max(abs(k),abs(sqrt(rho./G).*w)).^-1 ...
            .*exp(b/2*imag(kt)).*exp(b/2*imag(kl)) ...
            ; % antisymmetrical modes
end

% OUT-OF-PLANE MOTION
% quantities of interest
vinf = sqrt(xi).*vt ;
w0 = sqrt(12)/h*vinf ;
% wavenumbers (isotropic part mu_i² = k_i² + k²)
must2 = (w.^2-w0.^2)./vt.^2 ; % transverse shear
alpha = (vl./vinf).^2 + 1 ;
beta = 4.*(1-(w0./w).^2).*(vl./vinf).^2 ;
gamma = sqrt(alpha.^2 - beta) ;
mub2 = .5*(w./vl).^2.*(alpha + gamma) ; % bending
musa2 = .5*(w./vl).^2.*(alpha - gamma) ; % axial shear
% width-direction wavenumbers & trigo
kb = sqrt(mub2-k.^2) ; cb = cos(b./2.*kb) ; sb = sin(b./2.*kb) ; 
ksa = sqrt(musa2-k.^2) ; csa = cos(b./2.*ksa) ; ssa = sin(b./2.*ksa) ;
kst = sqrt(must2-k.^2) ; cst = cos(b./2.*kst) ; sst = sin(b./2.*kst) ;
% components of the system matrix
% Sst1 = -k ;
% Sst2 = kst.^2 - k.^2  ;
% Ast3 = (nu-1).*kst.*k ;
% Sb3 = w0.^2.*(nu.*k.^2 + kb.^2) ;
% Ssa3 = w0.^2.*(nu.*k.^2 + ksa.^2) ;
% Ab1 = kb.*(w.^2-vl.^2.*mub2) ;
% Ab2 = 2*kb.*k.*w0.^2 ;
% Asa1 = ksa.*(w.^2-vl.^2.*musa2) ;
% Asa2 = 2*ksa.*k.*w0.^2 ;
% precomputation
Dst = ...Ast3.*(Ab1.*Asa2 - Ab2.*Asa1) ...
      2*(1-nu).*kst.*kb.*ksa.*k.^2.*vl.^2.*(mub2-musa2) ...
      ; 
Dsa = ... Ssa3.*(Ab1.*Sst2-Ab2.*Sst1) ...
      kb.*(nu.*k.^2 + ksa.^2).*((w.^2-vl.^2.*mub2).*(kst.^2 - k.^2) + 2*k.^2.*w0.^2) ...
      ... w0.^2.*kb.*((nu-1).*k.^2 + musa2).*((w.^2-vl.^2.*mub2).*(must2 - 2*k.^2) + 2*k.^2.*w0.^2) ...
      ... w0.^2.*kb.*((nu-1).*k.^2 + musa2).*((w.^2-2*k.^2.*vt.^2).*must2 - vl.^2.*mub2.*(must2 - 2*k.^2)) ...
      ... w0.^2.*kb.*((nu-1).*k.^2 + musa2).*(w0.^2.*must2 + (vt.^2.*must2 - vl.^2.*mub2).*(must2 - 2*k.^2)) ...
      ; 
Db = ... Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
     ksa.*(nu.*k.^2 + kb.^2).*((w.^2-vl.^2.*musa2).*(kst.^2 - k.^2) + 2*k.^2.*w0.^2) ...
     ; 
% determinants
detBs = cst.*sb.*ssa.*Dst ...
      - csa.*sb.*sst.*Dsa ...
      + cb.*ssa.*sst.*Db ; % symmetric modes

detBa = sst.*cb.*csa.*Dst ...
      - ssa.*cb.*cst.*Dsa ...
      + sb.*csa.*cst.*Db ; % antisymmetric modes
% normalized
if normalize
    detBs = detBs...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^1.*ksa.^0.*gamma.^0)...
            ./abs(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % symmetrical modes
    detBa = detBa ...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^0.*ksa.^0.*gamma.^0)...
            ./(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % antisymmetrical modes
end

% complete determinant
D = detAs.*detAa.*detBa.*detBs ;

end
