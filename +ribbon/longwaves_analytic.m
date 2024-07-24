% %% LONG WAVES, IN PLANE MOTION
% clc
% clear all
% 
% syms V x B
% 
% %D = (1 + B*(V^2-1))*(V^2-2)^2 + 4*(x*V^2-1)*(1+B*(x*V^2-1)) ; % symmetric
% D = (1 + B*(x*V^2-1))*(V^2-2)^2 + 4*(V^2-1)*(1+B*(V^2-1)) ; % antisymmetric
% 
% D = expand(D/V^2) 
% 
% C = coeffs(D,V) ;
% a = C(3) ;
% b = C(2) ;
% c = C(1) ;
% % expand(a*V^4 + b*V^2 + c) - D
% 
% delta = b^2-4*a*c ;
% delta = expand(delta) ;
% Cdel = coeffs(delta,B)
% u = Cdel(2)/2
% v = expand(Cdel(3)-u^2)
% V = factor(v)
% 
% err = expand((1 + B*u)^2 + B^2*prod(V) - delta)
% 
% syms nu
% x = .5*(1-nu)
% deltanu = expand(subs(delta))
% simplify(sqrt(deltanu))
% 
% 
% %% LONG WAVES, OUT OF PLANE PLANE MOTION
% clc
% clear all
% 
% 
% % wavenumber expressions
% syms k V 
% syms x z s g w2 real positive
% kl2 = k^2*(x*V^2-1) ;
% kt2 = k^2*(V^2-1) ;
% kb2 = g.*k.^2.*V.^2 + z*k*V - k^2 ;
% ksa2 = g.*k.^2.*V.^2 - z*k*V - k^2 ;
% kst2 = kt2 - s^2 ;
% 
% % IN PLANE
% %D = 4*k^2*kl2*(1+w2*kl2) + (kt2-k^2)^2*(1+w2*kt2) ; % symmetric
% %D = 4*k^2*kt2*(1+w2*kt2) + (kt2-k^2)^2*(1+w2*kl2) ; % antisymmetric
% 
% 
% % OUT OF PLANE
% % WITH SIN/COSINES
% sincK = @(k2) ... 40*sin(.5*geo.b*sqrt(k2))./sqrt(k2) ...
%               (40 - 20*w2*k2 + 3*(w2*k2).^2) ...
%               ;
% cosK = @(k2) ... 8*cos(.5*geo.b*sqrt(k2)) ...
%               (8 - 12*w2*k2 + 3*(w2*k2).^2) ...
%             ;
% % D = 4.*x.*k.^2.*kb2.*ksa2.*(kb2-ksa2).*sincK(kb2).*sincK(ksa2).*cosK(kst2) ...
% %     + ksa2.*(kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*sincK(kst2).*sincK(ksa2).*cosK(kb2) ...
% %     - kb2.*(ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*sincK(kst2).*sincK(kb2).*cosK(ksa2) ...
% %     ; % symetric modes
% D = 4.*x.*k.^2.*kst2.*(kb2-ksa2).*sincK(kst2).*cosK(kb2).*cosK(ksa2) ...
%     + (kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*sincK(kb2).*cosK(ksa2).*cosK(kst2) ...
%     - (ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*sincK(ksa2).*cosK(kb2).*cosK(kst2) ...
%     ; % antisymetric modes
% 
% % Substitute
% D = subs(subs(D)) ;
% 
% % Remove common factors
% F = factor(D) ;
% D = F(end) 
% 
% 
% % Remove all powers of k>2
% [Ck,tk] = coeffs(D,k,'all') ;
% D = sum(Ck(end-2:end).*tk(end-2:end))
% 
% % Extract the coefficients with V
% [Cv,tv] = coeffs(D,V,'all') % [a*V^4 + b*V^2 + c]
% a = Cv(end-4)
% b = Cv(end-2)
% c = Cv(end)
% 
% % Determinant of the quadratic approximation
% %assumeAlso(w2*s^2>1) % needed for the square root simplification
% delta = b^2-4*a*c ;
% [Ckd,tkd] = coeffs(delta,k,'all') ;
% % approximate the square root of delta ~d0 + O(k)
% d0 = Ckd(end).*tkd(end) ; % zero order term
% %Ok = sum(Ckd(1:end-1).*tkd(1:end-1)) ; % O(k) term
% Ok = sum(Ckd(end-4:end-1).*tkd(end-4:end-1)) ; % O(k) term
% ratio = Ok./simplify(d0) ;
% sqrtdelta = sqrt(d0)*(1+ratio/2-ratio^2/8) 
% %[Cksd,tksd] = coeffs(sqrtdelta,k,'all') 
% %sqrtdelta = sum(Cksd(end-4:end).*tksd(end-4:end))
% 
% % Numerator (~=-b+sqrt(delta))
% numer = -b+sqrtdelta ;
% [Ckn,tkn] = coeffs(numer,k,'all') ;
% numer = sum(Ckn(end-4:end).*tkn(end-4:end)) ;
% 
% % Denominator of the solution (~=2a)
% twoA = 2*a ;
% [Ck2a,tk2a] = coeffs(twoA,k,'all') ;
% twoA = sum(Ck2a(end-2:end).*tk2a(end-2:end))
% 
% % Solutions
% V1 = simplify(expand(numer)/(twoA)) 
% %V2 = expand(-b-sqrtdelta)/(twoA)
% 
% syms h aa bb zeta nu real positive
% z = sqrt(12*x)/h ;
% s = sqrt(12*zeta)/h ;
% g = .5/zeta ;
% x = .5*(1-nu) ;
% w2 = aa^2*h^2/12 ; 
% % % D = subs(subs(subs(D))) 
% V1 = simplify(subs(subs(subs(V1))))
% 
% [cc,tt] = coeffs(V1,k)
% simplify(cc(end))


%% LONG WAVES, OUT OF PLANE PLANE MOTION
clc
clear all


% wavenumber expressions
syms V B % B = b²k²/12 = w2*k²
syms a real positive % a = b/h section aspect ratio
syms x real positive % x = vt²/vl² = (1-nu)/2
syms xi real positive % shear correction factor
syms g real positive % g = .5/xi
Bl = B*(x*V^2-1) ;
Bt = B*(V^2-1) ;
Bb = B.*V.^2/2/xi + a*sqrt(x)*sqrt(B)*V - B ;
Bsa = B.*V.^2/2/xi - a*sqrt(x)*sqrt(B)*V - B ;
Bst = Bt - a^2*xi ;

% SIN/COSINE
sincB = @(Bi) ... 40*sin(.5*geo.b*sqrt(k2))./sqrt(k2) ...
              (40 - 20*Bi + 0*3*(Bi).^2) ...
              ;
cosB = @(Bi) ... 8*cos(.5*geo.b*sqrt(k2)) ...
              (8 - 12*Bi + 0*3*(Bi).^2) ...
            ;

% IN PLANE
% D = 4*B*Bl*sincB(Bl)*cosB(Bt) + (Bt-B)^2*sincB(Bt)*cosB(Bl) ; % symmetric
% D = 4*B*Bt*sincB(Bt)*cosB(Bl) + (Bt-B)^2*sincB(Bl)*cosB(Bt) ; % antisymmetric

% OUT OF PLANE
% D = 4.*x.*B.*Bb.*Bsa.*(Bb-Bsa).*sincB(Bb).*sincB(Bsa).*cosB(Bst) ...
%     + Bsa.*(Bb+(1-2.*x).*B).*(2.*x.*B.*(Bt-Bst) + (Bl-Bsa).*(Bst-B)).*sincB(Bst).*sincB(Bsa).*cosB(Bb) ...
%     - Bb.*(Bsa+(1-2.*x).*B).*(2.*x.*B.*(Bt-Bst) + (Bl-Bb).*(Bst-B)).*sincB(Bst).*sincB(Bb).*cosB(Bsa) ...
%     ; % symetric modes
D = 4.*x.*B.*Bst.*(Bb-Bsa).*sincB(Bst).*cosB(Bb).*cosB(Bsa) ...
     + (Bb+(1-2.*x).*B).*(2.*x.*B.*(Bt-Bst) + (Bl-Bsa).*(Bst-B)).*sincB(Bb).*cosB(Bsa).*cosB(Bst) ...
     - (Bsa+(1-2.*x).*B).*(2.*x.*B.*(Bt-Bst) + (Bl-Bb).*(Bst-B)).*sincB(Bsa).*cosB(Bb).*cosB(Bst) ...
     ; % antisymetric modes

% Remove common factors
F = factor(D,V) ;
D = F(end) ;

% Remove all powers of B>1
[CB,tB] = coeffs(D,B,'all') ;
D = sum(CB(end-1:end).*tB(end-1:end)) ;

% Extract the coefficients with V
[Cv,tv] = coeffs(D,V,'all') ; % [p2*V^4 + p1*V^2 + p0]
p2 = Cv(end-4) ;
p1 = Cv(end-2) ;
p0 = Cv(end) ;

% Determinant of the quadratic approximation
delta = p1^2-4*p2*p0 ;
[CBd,tBd] = coeffs(delta,B,'all') ;
% approximate the square root of delta ~d0 + O(B)
d0 = simplify(CBd(end).*tBd(end)) ; % zero order term
OB = sum(CBd(end-2:end-1).*tBd(end-2:end-1)) ; % O(B) term
ratio = OB./d0 ; 
sqrtd0 = sqrt(d0) ;
sqrtdelta = sqrtd0*(1+ratio/2-ratio^2/8) ; 

% Numerator (~=-b(+/-)sqrt(delta))
[CBb,tBb] = coeffs(p1,B,'all') ;
if isequal(simplify(CBb(end)-sqrtd0),sym(0))
    numer = -p1+sqrtdelta ;
else
    numer = -p1-sqrtdelta ;
end
[CBn,tBn] = coeffs(numer,B,'all') ;
numer = sum(CBn(end-2:end).*tBn(end-2:end)) ;
numer = simplify(expand(numer)) ;

% Denominator of the solution (~=2a)
twoA = 2*p2 ;
[CB2a,tB2a] = coeffs(twoA,B,'all') ;
twoA = sum(CB2a(end-1:end).*tB2a(end-1:end)) ;

% Asymptotic Normalized velocity
Vb = numer/twoA ;

% Orders in B...
[cc,tt] = coeffs(Vb,B,'all') ;
cc = simplify(cc) ;

% Vb = V0 + B*V1
V0 = cc(end) % low-frequency asymptote
V1 = cc(end-1) % first-order term


%%
syms nu real 
simplify(subs(cc,x,(1-nu)/2))


