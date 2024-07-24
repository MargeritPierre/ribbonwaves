% Return all model matrices & elementary modes
[~,ik0] = min(abs(Kd(1,:))) ;
w = linspace(0,1,10000)'*2*pi*fc ; Wd(:,ik0) ; k = 0 ;

[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(geo,mat) ;
w2 = w.^2 ;
k2 = k.^2 ;

% IN-PLANE MOTION
% phase velocities
vl2 = Q./rho ;
vt2 = G./rho ;
% wavenumbers
kl = sqrt(w2./vl2 - k2) ; % longitudinal
kt = sqrt(w2./vt2 - k2) ; % transverse
% trigo
cl = cos(b./2.*kl) ; sl = sin(b./2.*kl) ; 
ct = cos(b./2.*kt) ; st = sin(b./2.*kt) ; 
% boundary conditions
Sl2 = (kt.^2-k2) ;
Al1 = 2*k.*kl ; 
St1 = (kt.^2-k2) ;
At2 = -2*k.*kt ; 
% determinants
Dips = Al1.*At2.*sl.*ct - Sl2.*St1.*st.*cl ; % symmetrical modes
Dipa = Al1.*At2.*st.*cl - Sl2.*St1.*sl.*ct ; % antisymmetrical modes
% normalized determinants
Dips_n = Dips./kt./(abs(w2)+eps) ; %(abs(kt)+eps) ; % symmetrical modes
Dipa_n = Dipa./kl./(abs(w2)+eps) ; %(abs(kl)+eps) ; % antisymmetrical modes

% OUT-OF-PLANE MOTION
% quantities of interest
vinf2 = xi.*vt2 ;
w02 = 12/h.^2.*vinf2 ;
% wavenumbers (isotropic part mu_i² = k_i² + k²)
must2 = (w2-w02)./vt2 ; % transverse shear
alpha = (vl2./vinf2) + 1 ;
beta2 = (1-(w02./w2)).*(vl2./vinf2) ;
gamma = sqrt(alpha.^2 - 4*beta2) ;
mub2 = .5*(w2./vl2).*(alpha + gamma) ; % bending
musa2 = .5*(w2./vl2).*(alpha - gamma) ; % axial shear
% width-direction wavenumbers
kb = sqrt(mub2-k2) ;
ksa = sqrt(musa2-k2) ;
kst = sqrt(must2-k2) ;
% trigo
cb = cos(b./2.*kb) ; sb = sin(b./2.*kb) ; 
csa = cos(b./2.*ksa) ; ssa = sin(b./2.*ksa) ;
cst = cos(b./2.*kst) ; sst = sin(b./2.*kst) ;
% boundary conditions
Sst1 = -k ;
Sst2 = kst.^2 - k2  ;
Ast3 = (nu-1).*kst.*k ;
Sb3 = w02.*(nu.*k2 + kb.^2) ;
Ssa3 = w02.*(nu.*k2 + ksa.^2) ;
Ab1 = kb.*(w2-vl2.*mub2) ;
Ab2 = 2*kb.*k.*w02 ;
Asa1 = ksa.*(w2-vl2.*musa2) ;
Asa2 = 2*ksa.*k.*w02 ;
% precomputation
Dst = 2*kst.*kb.*ksa.*k2.*(mub2-musa2); 
Dsa = kb.*(ksa.^2 + nu.*k2).*(k2.*(kt.^2-kst.^2) + (1./(1-nu)).*(kl.^2-kb.^2).*(kst.^2-k2) ) ; 
Db =  ksa.*(kb.^2 + nu.*k2).*(k2.*(kt.^2-kst.^2) + (1./(1-nu)).*(kl.^2-ksa.^2).*(kst.^2-k2) ) ; 
% determinants
Dops = cst.*sb.*ssa.*Dst ...
      - csa.*sb.*sst.*Dsa ...
      + cb.*ssa.*sst.*Db ; % symmetric modes
Dopa = sst.*cb.*csa.*Dst ...
      - ssa.*cb.*cst.*Dsa ...
      + sb.*csa.*cst.*Db ; % antisymmetric modes
% normalized determinants
Dops_n = Dops./must2./kst./(abs(w2)+eps) ; %(abs(must2)+eps) ; % symmetrical modes
Dopa_n = Dopa./must2./kb./ksa./(abs(w2)+eps) ; %(abs(must2)+eps)./(abs(kb)+eps)./(abs(ksa)+eps) ; % antisymmetrical modes

BpA = .5*b.*w./sqrt(vl2).*sqrt(alpha+2.*sqrt(beta2)) ;
BmA = .5*b.*w./sqrt(vl2).*sqrt(alpha-2.*sqrt(beta2)) ;
pp = (1-sqrt(beta2)).*sqrt(alpha+2.*sqrt(beta2)) ;
mm = (1+sqrt(beta2)).*sqrt(alpha-2.*sqrt(beta2)) ;

Dops_a = sin(BpA).*mm - sin(BmA).*pp ;
Dopa_a = sin(BpA).*mm + sin(BmA).*pp ;

clf 
plot(w./sqrt(w02),log10(abs(Dops_n./sst))+14,'r') ; 
plot(w./sqrt(w02),log10(abs(Dops_a)),':r')
plot(w./sqrt(w02),log10(abs(Dopa_n./cst))+11,'b') ; 
plot(w./sqrt(w02),log10(abs(Dopa_a))-3,':b')

%% Small ratios w/w0
A = sqrt(sqrt(3)/2)*b.*sqrt(w./sqrt(vl2)/h) ;
Dops_aa = tan(A) + tanh(A) ;
Dopa_aa = tan(A) - tanh(A) ;
plot(w./sqrt(w02),log10(abs(Dops_aa)),'--r') ;
plot(w./sqrt(w02),log10(abs(Dopa_aa))-3,'--b') ;


%% Large ratios w/w0
Dops_aa = sin(.5*b.*w.*(1./sqrt(vl2)+1./sqrt(vinf2))) - sin(.5*b.*w.*(1./sqrt(vl2)-1./sqrt(vinf2))) ;
Dopa_aa = sin(.5*b.*w.*(1./sqrt(vl2)+1./sqrt(vinf2))) + sin(.5*b.*w.*(1./sqrt(vl2)-1./sqrt(vinf2))) ;
plot(w./sqrt(w02),log10(abs(Dops_aa)),'-.r')
plot(w./sqrt(w02),log10(abs(Dopa_aa))-3,'-.b')
n = 1:100 ;
wn_inf = n*pi/b*sqrt(vinf2)
wn_l = n*pi/b*sqrt(vl2)
% plot(wn_inf./sqrt(w02).*[1;1],get(gca,'ylim'),'k')
% plot(wn_l./sqrt(w02).*[1;1],get(gca,'ylim'),'m')


