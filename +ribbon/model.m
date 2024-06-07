function m = model(GEO,MAT,k,w)
% Return all model matrices & elementary modes
[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(GEO,MAT) ;
w2 = w.^2 ;
k2 = k.^2 ;

% IN-PLANE MOTION
% phase velocities
vl2 = Q./rho ;
vt2 = G./rho ;
% wavenumbers
kl = sqrt(w2./vl2 - k2) ; % longitudinal
kt = sqrt(w2./vt2 - k2) ; % transverse
% polarizations
Ulp = {k ; kl} ;
Ulm = {k ; -kl} ;
Utp = {kt ; -k} ;
Utm = {-kt ; -k} ;
% trigo
cl = cos(b./2.*kl) ; sl = sin(b./2.*kl) ; 
ct = cos(b./2.*kt) ; st = sin(b./2.*kt) ; 
% boundary conditions
Sl2 = (kt.^2-k2) ;
Al1 = 2*k.*kl ; 
St1 = (kt.^2-k2) ;
At2 = -2*k.*kt ; 
% Matrices
Mips = {sl.*Al1 st.*St1 ; ...
        cl.*Sl2 ct.*At2} ;
Mipa = {cl.*Al1 ct.*St1 ; ...
        sl.*Sl2 st.*At2} ;
% determinants
Dips = Al1.*At2.*sl.*ct - Sl2.*St1.*st.*cl ; % symmetrical modes
Dipa = Al1.*At2.*st.*cl - Sl2.*St1.*sl.*ct ; % antisymmetrical modes
% normalized determinants
Dips_n = Dips./kt./(abs(w2)+eps) ; %(abs(kt)+eps) ; % symmetrical modes
Dipa_n = Dipa./kl./(abs(w2)+eps) ; %(abs(kl)+eps) ; % antisymmetrical modes
% amplitudes
[als,ata] = eigenvectors(Mips{:}) ;
[ala,ats] = eigenvectors(Mipa{:}) ;

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
% polarizations
Ust = {0*k ; kst ; -k} ;
Ub = {vl2.*mub2 + w02 - w2 ; 1i*w02.*k ; 1i.*w02.*kb} ;
Usa = {vl2.*musa2 + w02 - w2 ; 1i*w02.*k ; 1i.*w02.*ksa} ;
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
% matrices
Mops = {sst.*Sst1 sb.*Ab1 ssa.*Asa1 ; ...
        sst.*Sst2 sb.*Ab2 ssa.*Asa2 ; ...
        cst.*Ast3 cb.*Sb3 csa.*Ssa3} ;
Mopa = {cst.*Sst1 cb.*Ab1 csa.*Asa1 ; ...
        cst.*Sst2 cb.*Ab2 csa.*Asa2 ; ...
        sst.*Ast3 sb.*Sb3 ssa.*Ssa3} ;
% precomputation
Dst = ... Ast3.*(Ab1.*Asa2 - Ab2.*Asa1) ...
      ...2*(1-nu).*kst.*kb.*ksa.*k2.*vl2.*(mub2-musa2) ...
      4*kst.*kb.*ksa.*k2.*vt2.*vl2.*(mub2-musa2) ...
      ; 
Dsa = ... Ssa3.*(Ab1.*Sst2-Ab2.*Sst1) ...
      ...kb.*(nu.*k2 + ksa.^2).*((w2-vl2.*mub2).*(kst.^2 - k2) + 2*k2.*w02) ...
      ...kb.*(vl2.*musa2 - 2*vt2.*k2).*( must2.*(w2 - vl2.*mub2) + 2*k2.*(vl2.*mub2 - vt2.*must2) ) ...
      ...kb.*(vl2.^2.*mub2.*musa2.*(2*k2-must2) + 2*vl2.*vt2.*k2.*must2.*(mub2+musa2)+2*vt2.*k2.*must2.*(2*vt2.*k2-w2)-4.*vl2.*vt2.*k2.^2.*mub2 + vl2.*must2.*musa2.*(w2-4*vt2.*k2)) ...
      kb.*(vl2.*w2.*must2.*(4*k2-must2)/xi + 4*vt2.^2.*k2.^2.*must2 -4.*vl2.*vt2.*k2.^2.*mub2 + vl2.*must2.*musa2.*(w2-4*vt2.*k2)) ...
      ; 
Db =  ... Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
      ...ksa.*(nu.*k2 + kb.^2).*((w2-vl2.*musa2).*(kst.^2 - k2) + 2*k2.*w02) ...
      ...ksa.*(vl2.*mub2 - 2*vt2.*k2).*( must2.*(w2 - vl2.*musa2) + 2*k2.*(vl2.*musa2 - vt2.*must2) ) ...
      ...ksa.*(vl2.^2.*mub2.*musa2.*(2*k2-must2) + 2*vl2.*vt2.*k2.*must2.*(mub2+musa2)+2*vt2.*k2.*must2.*(2*vt2.*k2-w2)-4.*vl2.*vt2.*k2.^2.*musa2 + vl2.*must2.*mub2.*(w2-4*vt2.*k2)) ...
      ksa.*(vl2.*w2.*must2.*(4*k2-must2)/xi + 4*vt2.^2.*k2.^2.*must2 -4.*vl2.*vt2.*k2.^2.*musa2 + vl2.*must2.*mub2.*(w2-4*vt2.*k2)) ...
     ; 
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
  
% Build the model structure
m = [] ;
for field = [...
                "vl2","vt2","kl","kt" ...
                ...,"Ul","Ut","cl","sl","ct","st","Sl2","Al1","St1","At2" ...
                ,"Dips","Dipa","Dops","Dopa" ...
                ,"Dips_n","Dipa_n","Dops_n","Dopa_n" ...
                ...,"als","ata","ala","ats" ...
                ]
            m.(field) = eval(field) ;
end

end

function varargout = eigenvectors(varargin)
% Compute eigenvectors v s.t A.v = 0
v = repmat({NaN(size(varargin{1}))},[sqrt(numel(varargin)) 1]) ;
switch nargin
    case 4 % 2x2 matrix [a c ; b d]
        [a,b,c,d] = deal(varargin{:}) ;
    % Norms of the two vectors 
        n1 = sqrt(abs(a).^2 + abs(c).^2) ;
        n2 = sqrt(abs(b).^2 + abs(d).^2) ;
        geq = n1>=n2 ; ll = ~geq ;
        v{1}(geq) = c(geq)./n1(geq) ;
        v{2}(geq) = -a(geq)./n1(geq) ;
        v{1}(ll) = d(ll)./n2(ll) ;
        v{2}(ll) = -b(ll)./n2(ll) ;
    case 9 % 3x3 matrix
end
varargout = v ;
end

