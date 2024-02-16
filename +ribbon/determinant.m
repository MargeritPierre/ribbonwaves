%% DETERMINANT COMPUTATION
function [D,detAs,detAa,detBs,detBa] = determinant(h,b,Q,G,rho,w,k,normalize,xi)
if nargin<8 || isempty(normalize) ; normalize = false ; end
nu = 1-2.*G./Q ;

% IN-PLANE MOTION
% wavenumbers
kl = sqrt(rho.*w.^2./Q - k.^2) ; % longitudinal
kt = sqrt(rho.*w.^2./G - k.^2) ; % transverse
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
            .*exp(b/2*imag(kt+kl)) ; % symmetrical modes
    detAa = detAa./Al...
            .*abs(k).^2.*abs(w).^-2 ...
            .*abs(kt).^1.*max(abs(k),abs(sqrt(rho./G).*w)).^-1 ...
            ....*kl.^-1.*max(k,sqrt(rho/Q).*w).^1 ...
            .*exp(b/2*imag(kt)).*exp(b/2*imag(kl)) ...
            ; % antisymmetrical modes
end

% OUT-OF-PLANE MOTION
if nargin<9 ; xi = sqrt(pi^2/12) ; end % shear correction factor
eta = 12./(h.^2) ; % shape factor
% wavenumbers (isotropic part mui2 = ki2 + k2)
ssss = sqrt((1./xi./G - 1./Q).^2 + 4*eta./rho./Q./(w.^2)) ; 
mub2 = rho.*w.^2./2.*((1./xi./G + 1./Q) + ssss) ; % bending
musa2 = rho.*w.^2./2.*((1./xi./G + 1./Q) - ssss) ; % axial shear
must2 = rho.*w.^2./G-eta.*xi ; % transverse shear
% width-direction wavenumbers & trigo
kb = sqrt(mub2-k.^2) ; cb = cos(b./2.*kb) ; sb = sin(b./2.*kb) ; 
ksa = sqrt(musa2-k.^2) ; csa = cos(b./2.*ksa) ; ssa = sin(b./2.*ksa) ;
kst = sqrt(must2-k.^2) ; cst = cos(b./2.*kst) ; sst = sin(b./2.*kst) ;
% components of the system matrix
% Sb3 = eta.*xi.*G.*((nu-1)*k.^2 + mub2) ;
% Ssa3 = eta.*xi.*G.*((nu-1)*k.^2 + musa2) ;
% Ab1 = kb.*(rho.*w.^2-Q.*mub2) ;
% Ab2 = 2.*eta.*xi.*G.*k.*kb ;
% Asa1 = ksa.*(rho.*w.^2-Q.*musa2) ;
% Asa2 = 2.*eta.*xi.*G.*k.*ksa ;
% Sst1 = k ;
% Sst2 = k.^2 - kst.^2 ;
% Ast3 = (1-nu).*k.*kst ;
% precomputation
Ast3Ab1Asa2mAst3Ab2Asa1 = ...Ast3.*(Ab1.*Asa2-Ab2.*Asa1) ...
                          -4.*G.*rho.*w.^2.*k.^2.*kst.*kb.*ksa.*ssss ...
                          ; 
Ssa3Ab2Sst1mSsa3Ab1Sst2 = ...Ssa3.*(Ab2.*Sst1-Ab1.*Sst2) ...
                          kb.*(ksa.^2 + nu.*k.^2).*((rho.*w.^2-2.*G.*k.^2).*(k.^2+kst.^2) + Q.*(k.^2 + kb.^2).*(k.^2-kst.^2)) ...
                          ; 
Sb3Asa1Sst2mSb3Asa2Sst1 = ...Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
                          -ksa.*(kb.^2 + nu.*k.^2).*((rho.*w.^2-2.*G.*k.^2).*(k.^2+kst.^2) + Q.*(k.^2 + ksa.^2).*(k.^2-kst.^2)) ...
                          ; 
% determinants
detBs = cst.*sb.*ssa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + csa.*sb.*sst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + cb.*ssa.*sst.*Sb3Asa1Sst2mSb3Asa2Sst1 ; % symmetric modes

detBa = sst.*cb.*csa.*Ast3Ab1Asa2mAst3Ab2Asa1 ...
      + ssa.*cb.*cst.*Ssa3Ab2Sst1mSsa3Ab1Sst2 ...
      + sb.*csa.*cst.*Sb3Asa1Sst2mSb3Asa2Sst1 ; % antisymmetric modes
% normalized
if normalize
    detBs = detBs...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^1.*ksa.^0.*ssss.^0)...
            ./abs(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % symmetrical modes
    detBa = detBa ...
            ./abs(w.^3.*k.^0.*kst.^0.*kb.^0.*ksa.^0.*ssss.^0)...
            ./(k.^2+kst.^2) ...
            .*exp(b/2*(imag(kb+kst)-abs(imag(ksa)))) ...
            ; % antisymmetrical modes
end

% complete determinant
D = detAs.*detAa.*detBa.*detBs ;

end
