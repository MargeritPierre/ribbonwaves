% Return all model matrices & elementary modes
kr = .2*linspace(-1,1,500)*kc ;
ki = 2*linspace(-1,1,500)*kc ;
[KR,KI] = meshgrid(kr,ki) ;
k = KR + 1i*KI ;
w = 2*pi*fc*1e-7 ;

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


% FIRST-ORDER APPROXIMATION
THETA = (12*xi/h^2)./k2 ;
DELTA = sqrt(3/h^2).*(w./k2./sqrt(vl2)) ;

kst_a = 1i*k.*sqrt(1+THETA) ;
kb_a = 1i*k.*(1-DELTA) ;
ksa_a = 1i*k.*(1+DELTA) ;

Dst_a = -8i*k.^7.*DELTA.*sqrt(1+THETA) ;
Dsa_a = 2i*k.^7.*(-DELTA.*(2*THETA+2-.5*(1-nu).*THETA)-.5*(1-nu).*THETA) ;
Db_a = 2i*k.^7.*(DELTA.*(2*THETA+2-.5*(1-nu).*THETA)-.5*(1-nu).*THETA) ;
DsamDb_a = -8i.*k.^7.*(DELTA.*(THETA+1-.25*(1-nu).*THETA)) ;
DsapDb_a = -8i.*k.^7.*(.25*(1-nu).*THETA) ;
sst_a = 1i*sinh(0.5*b*k.*sqrt(1+THETA)) ;
cst_a = cosh(0.5*b*k.*sqrt(1+THETA)) ;
sb_a = 1i*(sinh(.5*b*k)-.5*b.*k.*DELTA.*cosh(.5*b*k)) ; 
cb_a = cosh(.5*b*k)-.5*b.*k.*DELTA.*sinh(.5*b*k) ;
ssa_a = 1i*(sinh(.5*b*k)+.5*b.*k.*DELTA.*cosh(.5*b*k)) ; 
csa_a = cosh(.5*b*k)+.5*b.*k.*DELTA.*sinh(.5*b*k) ;
sbssa_a = .5*(1-cosh(b.*k)) ;
cbcsa_a = .5*(1+cosh(b.*k)) ;
sbcsa_a = .5i*(sinh(b*k)-b*k.*DELTA) ;
cbssa_a = .5i*(sinh(b*k)+b*k.*DELTA) ;

Dops_a = cosh(0.5*b.*k.*sqrt(1+THETA)).*(1-cosh(b.*k)).*sqrt(1+THETA) ...
         + sinh(0.5*b.*k.*sqrt(1+THETA)).*(...
                sinh(b*k).*(THETA+1)...
                - .25*(1-nu).*THETA.*(sinh(b*k)+b.*k)...
            ) ...
         ...cst_a.*sbssa_a.*Dst_a - sbcsa_a.*sst_a.*Dsa_a + cbssa_a.*sst_a.*Db_a ...
        ; % symmetric modes

Dopa_a = sinh(0.5*b.*k.*sqrt(1+THETA)).*(1+cosh(b.*k)).*sqrt(1+THETA) ...
            - cosh(0.5*b*k.*sqrt(1+THETA)).*(...
                sinh(b*k).*(THETA+1) ...
                - .25*(1-nu).*THETA.*(sinh(b*k)-(b.*k)) ...
            ) ...
        ... sst_a.*cbcsa_a.*Dst_a - cbssa_a.*cst_a.*Dsa_a + sbcsa_a.*cst_a.*Db_a ...
        ; % symmetric modes
    
a = b/h ;
phi = a*sqrt(3*xi) ;
y = k.*b./(2*phi) ;

Dops_a_norm = cosh(phi.*sqrt(1+y.^2)).*(1-cosh(2*phi.*y)).*y.*sqrt(1+y.^2) ...
         + sinh(phi.*sqrt(1+y.^2)).*(...
                sinh(2*phi.*y).*(y.^2+1)...
                - .25*(1-nu).*(sinh(2*phi.*y)+2*phi.*y)...
            ) ...
         ...cst_a.*sbssa_a.*Dst_a - sbcsa_a.*sst_a.*Dsa_a + cbssa_a.*sst_a.*Db_a ...
        ; % symmetric modes

Dopa_a_norm = sinh(phi.*sqrt(1+y.^2)).*(1+cosh(2*phi.*y)).*y.*sqrt(1+y.^2) ...
            - cosh(phi.*sqrt(1+y.^2)).*(...
                sinh(2*phi.*y).*(y.^2+1) ...
                - .25*(1-nu).*(sinh(2*phi.*y)-(2*phi.*y)) ...
            ) ...
        ... sst_a.*cbcsa_a.*Dst_a - cbssa_a.*cst_a.*Dsa_a + sbcsa_a.*cst_a.*Db_a ...
        ; % symmetric modes
  
I = kst ; I_a = kst_a ; 
I = kb ; I_a = kb_a ;  
I = ksa ; I_a = ksa_a ;  
I = Dst ; I_a = Dst_a ; 
I = Db ; I_a = Db_a ; 
I = Dsa ; I_a = Dsa_a ; 
I = cb.*ssa ; I_a = cbssa_a ; 
I = Dops ; I_a = -4i.*k.^7.*DELTA.*Dops_a ; 
I = Dopa ; I_a = 4.*k.^7.*DELTA.*Dopa_a ; 
I = Dops ; I_a = 4.*k.^7.*DELTA.*Dops_a_norm./y.^2 ; 
I = Dopa ; I_a = 4.*k.^7.*DELTA.*Dopa_a_norm./y.^2 ; 
%I = THETA ; I_a = (1./x.^2) ; 

clf reset
for ii = {I I_a (I_a-I)./abs(I) (abs(I_a)-abs(I))./abs(I)}
mysubplot(2,2,numel(get(gcf,'children'))/2+1) ;
    iii = ii{end} ;
    %iii = iii./medfilt2(abs(iii),[3 3]) ;
    surf(KR/kc,KI/kc,log10(abs(iii))) ;
    shading interp
    %colormap gray
    %caxis([-.3 0])
    axis tight
    colorbar%('handlevisibility','off') ;
end

%% ROOT FINDING
% We search the values y as the roots of Dop(s/a)_a_norm
% We approximate the different terms with their respective series
% expansion: https://en.wikipedia.org/wiki/Taylor_series
% fi(x) ~ Pi(x) = sum_n^N{pi.*x.^n}
N = 150 ; % polynomial order approximation
n = 0:N-1 ;
% SQRT(1+y²): see https://planetmath.org/taylorexpansionofsqrt1x
P_sqrt_Y2p1 = Pbuild(2*n,[1 .5 (-1).^(n(3:end)-1) .*factorial(2*n(3:end)-3)./4.^(n(3:end)-1)./factorial(n(3:end))./factorial(n(3:end)-2)],N) ;
% SINH(x),COSH(x) ;
% SINH(x) = x.^(2*n+1)./factorial(2*n+1)
P_sinh_X = Pbuild(2*n+1,1./factorial(2*n+1),N) ;
% COSH(x) = x.^(2*n)./factorial(2*n)
P_cosh_X = Pbuild(2*n,1./factorial(2*n),N) ;
% SINH(2*phi*y)
P_sinh_2phiY = Pcomp(P_sinh_X,2*phi*[0;1],N) ;
% COSH(2*phi*y)
P_cosh_2phiY = Pcomp(P_cosh_X,2*phi*[0;1],N) ;
% SINH(phi*sqrt(1+y²)) = SINH(x) with x = phi*sqrt(1+y²)
P_sinh_phisqrt_Y2p1 = Pcomp(P_sinh_X,phi*P_sqrt_Y2p1,N) ;
% COSH(phi*sqrt(1+y²)) = COSH(x) with x = phi*sqrt(1+y²)
P_cosh_phisqrt_Y2p1 = Pcomp(P_cosh_X,phi*P_sqrt_Y2p1,N) ;
% Dops_a_norm = cosh(phi.*sqrt(1+y.^2)).*(1-cosh(2*phi.*y)).*y.*sqrt(1+y.^2)
P_ops_a = Ptimes(Ptimes(Ptimes(P_cosh_phisqrt_Y2p1,Pplus(1,-P_cosh_2phiY)),[0;1]),P_sqrt_Y2p1,N) ...
...+ sinh(phi.*sqrt(1+y.^2)).*( 
    + Ptimes( P_sinh_phisqrt_Y2p1, ...
...   sinh(2*phi.*y).*(y.^2+1)
        Ptimes(P_sinh_2phiY,[1 0 1],N) ...
...   - .25*(1-nu).*(sinh(2*phi.*y)+2*phi.*y)
        -.25.*(1-nu).*Pplus(P_sinh_2phiY,[0 2*phi]) ...
... ) 
    ,N) ...
    ;

% Display tests
ytest = 1*linspace(-1,1,1000) ; 
Ytest = ytest(:).^n ;
Dops_test =  cosh(phi.*sqrt(1+ytest.^2)).*(1-cosh(2*phi.*ytest)).*ytest.*sqrt(1+ytest.^2) ...
         + sinh(phi.*sqrt(1+ytest.^2)).*(...
                sinh(2*phi.*ytest).*(ytest.^2+1)...
                - .25*(1-nu).*(sinh(2*phi.*ytest)+2*phi.*ytest)...
            ) ...
         ...cst_a.*sbssa_a.*Dst_a - sbcsa_a.*sst_a.*Dsa_a + cbssa_a.*sst_a.*Db_a ...
        ; % symmetric modes
clf ; plot(ytest,sqrt(1+ytest.^2)) ; plot(ytest,Ytest*P_sqrt_Y2p1(:),':') ;
clf ; plot(ytest,sinh(2*phi.*ytest)) ; plot(ytest,Ytest*P_sinh_2phiY(:),':') ;
clf ; plot(ytest,cosh(2*phi.*ytest)) ; plot(ytest,Ytest*P_cosh_2phiY(:),':') ;
clf ; plot(ytest,sinh(phi*sqrt(1+ytest.^2))) ; plot(ytest,Ytest*P_sinh_phisqrt_Y2p1(:),':') ;
clf ; plot(ytest,cosh(phi*sqrt(1+ytest.^2))) ; plot(ytest,Ytest*P_cosh_phisqrt_Y2p1(:),':') ;
clf ; plot(ytest,Dops_test) ; plot(ytest,Ytest*P_ops_a(:),':') ;

% PLOT ROOTS ON THE DETERMINANT
yr_ops = roots(flip(P_ops_a)) ;
kr_ops = 2*phi*yr_ops/b ;
clf
iii = Dops_a_norm ;
%iii = iii./medfilt2(abs(iii),[3 3]) ;
srf = surf(KR/kc,KI/kc,log10(abs(iii))) ;
srf.CData = srf.CData ; srf.ZData = srf.ZData*0 ;
shading interp
%colormap gray
%caxis([-.3 0])
axis tight
colorbar%('handlevisibility','off') ;
plot(real(kr_ops)/kc,imag(kr_ops)/kc,'o') ;

%% UTIL FUNCTIONS
function P = Pbuild(n,an,N)
% Build a polynomial 
    P = full(sparse(1,n+1,an)) ;
    if nargin>=3 ; P = P(1:N) ; end
end
function Pc = Ptimes(Pa,Pb,N)
% Product of two polynomials
    Pc = conv(Pa,Pb) ;
    if nargin>=3 ; Pc = Pc(1:N) ; end
end

function Pc = Pplus(Pa,Pb)
% Sum of two polynomials
    n = max(numel(Pa),numel(Pb)) ;
    Pa(end+1:n) = 0 ; Pb(end+1:n) = 0 ; 
    Pc = Pa+Pb ;
end

function Pc = Pcomp(Pa,Pb,N)
% Polynomial composition so that C(x) = A(B(x))
    Na = numel(Pa) ; Nb = numel(Pb) ;
    if nargin<3 ; N = (Na-1)*Nb ; end
    Pbn = [1 zeros(1,N-1)] ; % unitary polynomial
    Pc = zeros(1,N) ;
    for nn = 1:Na
        Pc(1:numel(Pbn)) = Pc(1:numel(Pbn)) + Pa(nn).*Pbn ;
        Pbn = Ptimes(Pbn,Pb,N) ; % n-th power of B
    end
end

