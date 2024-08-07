%% WAVENUMBER SPECTRUM OF A PRISMATIC WAVEGUIDE
% Comparison between the SAFE method and the reduced ribbon model
clc
clear all
% Computed variable
    toCompute = 'k' ; % 'k' or 'w'
% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*10 ; % width (mm)
    dx = min([geo.b geo.h])/5 ; % SAFE element size
% Material
    mat = [] ;
    mat.E = 70e3*(1+0.0i) ; % Young modulus (MPa)
    mat.nu = 30/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
% Adimensionalized numbers
    fc = .5*sqrt(real(mat.G)/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Wavenumber & frequency domains
    kp = kc*[-.25i 0 1.1] ; % wavenumber keypoints
    wp = [0 1.1]*2*pi*fc ; % frequency range
% Computation points
    nPts = 100 ;
    switch toCompute % re-interpolate from keypoints
        case 'w'
            kr = interp1([0 cumsum(abs(diff(kp)))],kp,linspace(0,1,nPts)*sum(abs(diff(kp)))) ;
        case 'k'
            wr = interp1([0 cumsum(abs(diff(wp)))],wp,linspace(0,1,nPts)*sum(abs(diff(wp)))) ;
    end
% RIBBON MODEL
    % Map of the determinant
        nPix = 500 ;
        [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,kp,wp,nPix) ;
    % Direct computations
        switch toCompute
            case 'w'
                [Wips,Wipa,Wops,Wopa] = ribbon.computeW(geo,mat,kr) ;
                Kr = repmat(kr,[size(Wips,1) 1]) ;
            case 'k'
                [Kips,Kipa,Kops,Kopa] = ribbon.computeK(geo,mat,wr) ;
                Wr = repmat(wr,[size(Kips,1) 1]) ;
        end
% SAFE
    % Eigenvalue options
        nModes = 25 ;
    % Build the mesh
        mesh = safe.mesh.gridmesh([geo.b geo.h],dx) ; % quad mesh
        mesh = safe.mesh.quad8(mesh) ; % serendipity element
        mesh = safe.mesh.quadrature(mesh,'quad2') ; % generate a quadrature rule
        clf ; axis equal ; safe.mesh.plotmesh(mesh) ; 
    % Compute the SAFE frequencies & modes
    switch toCompute
        case 'w'
            [Us,Ws] = safe.computeW(mesh,mat,kr,nModes,false) ;
            Ks = repmat(kr,[nModes 1]) ;
        case 'k'
            [Us,Ks] = safe.computeK(mesh,mat,wr,nModes,false) ;
            Ws = repmat(wr,[nModes 1]) ;
    end
%% Display
    clf ; hold on
    % For legend only
        plot(NaN,NaN,'color',colors(1,:),'displayname','IP/S') ;
        plot(NaN,NaN,'color',colors(2,:),'displayname','IP/A') ;
        plot(NaN,NaN,'color',colors(3,:),'displayname','OP/S') ;
        plot(NaN,NaN,'color',colors(4,:),'displayname','OP/A') ;
    % Det map
        surf(real(Kd)+imag(Kd),Wd,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    switch toCompute
        case 'w'
        % Ribbon model solution
            plot3(real(Kr(:))+imag(Kr(:)),real(Wips(:)),imag(Wips(:)),'.','color',colors(1,:),'displayname','$\omega_{IP}^S$') ; 
            plot3(real(Kr(:))+imag(Kr(:)),real(Wipa(:)),imag(Wipa(:)),'.','color',colors(2,:),'displayname','$\omega_{IP}^A$') ; 
            plot3(real(Kr(:))+imag(Kr(:)),real(Wops(:)),imag(Wops(:)),'.','color',colors(3,:),'displayname','$\omega_{OP}^S$') ; 
            plot3(real(Kr(:))+imag(Kr(:)),real(Wopa(:)),imag(Wopa(:)),'.','color',colors(4,:),'displayname','$\omega_{OP}^A$') ; 
        % SAFE solutions
            plot3(real(Ks(:))+imag(Ks(:)),real(Ws(:)),imag(Ws(:)),'ok','displayname','SAFE','markersize',4,'linewidth',.5) ;  
        case 'k'
        % Ribbon model solution
            plot3(real(Kips(:)),real(Wr(:))+imag(Wr(:)),imag(Kips(:)),'.','color',colors(1,:),'displayname','$\omega_{IP}^S$') ; 
            plot3(real(Kipa(:)),real(Wr(:))+imag(Wr(:)),imag(Kipa(:)),'.','color',colors(2,:),'displayname','$\omega_{IP}^A$') ; 
            plot3(real(Kops(:)),real(Wr(:))+imag(Wr(:)),imag(Kops(:)),'.','color',colors(3,:),'displayname','$\omega_{OP}^S$') ; 
            plot3(real(Kopa(:)),real(Wr(:))+imag(Wr(:)),imag(Kopa(:)),'.','color',colors(4,:),'displayname','$\omega_{OP}^A$') ; 
        % SAFE solutions
            plot3(real(Ks(:)),real(Ws(:))+imag(Ws(:)),imag(Ks(:)),'ok','displayname','SAFE','markersize',4,'linewidth',.5) ;  
    end
% Legends
    xlabel 'Wavenumber $k$ (rad/mm)' ; 
    ylabel 'Frequency $\omega$ (rad/s)'
    lgd = legend ; lgd.Location = 'southeast' ;
% Rescaling
    obj = allchild(gca) ;
    [xdata,ydata,zdata] = deal(get(obj,'XData'),get(obj,'YData'),get(obj,'ZData')) ;
    set(obj,{'xdata'},cellfun(@(x)x./kc,xdata,'uni',false))
    set(obj,{'ydata'},cellfun(@(y)y./2/pi/fc,ydata,'uni',false))
    switch toCompute
        case 'w'
            set(obj,{'zdata'},cellfun(@(z)z./2/pi/fc,zdata,'uni',false))
        case 'k'
            set(obj,{'zdata'},cellfun(@(z)z./kc,zdata,'uni',false))
    end
    xlabel 'Normalized Wavenumber $\frac{k \times h}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2h}{c_t}$'
    set(gca,'xlim',[min(real(kp)+imag(kp)) max(real(kp)+imag(kp))]./kc,'ylim',[min(wp) max(wp)]./2./pi./fc) ;

    
%% CUTOFF FREQS
[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(geo,mat) ;
wt = linspace(0,5,10000)*2*pi*fc ;
m = ribbon.model(geo,mat,0,wt) ;
x = G./Q ;
a = b./h ;
alpha = 1+(1./xi./x) ;
beta = sqrt((alpha-1).*(1-(12./h^2)*xi.*(G./rho)./wt.^2)) ;
vl = sqrt(Q./rho) ;

m = -1 ; % symmetry

sinc = @(x)sin(x)./x ;
d0 = sinc(.5*b*wt./vl.*sqrt(alpha+2*beta)).*(1+beta) - m*sinc(.5*b*wt./vl.*sqrt(alpha-2*beta)).*(1-beta) ;

sinc_sqrt = @(X)math.poly.sincX(X.N).keeporders(2).comp(X,X.N) ;
switch 3
    case 0 % very low frequency approximation
        n = 1:3 ;
        Wr = pi^2/sqrt(3)*(n-m/4).^2/a/sqrt(x) ;
    case 1 % low frequency approximation
        W = math.poly.X(50) ;
        W2 = W*W ;
        WBeta = 1i*sqrt(3.*xi.*a.^2-W2)./sqrt(xi.*x) ;
        Dw = sinc_sqrt(x*(W2*alpha+2*W*WBeta)).*(W+WBeta) - m*sinc_sqrt(x*(W2*alpha-2*W*WBeta)).*(W-WBeta) ;
        Dw.coeffs = Dw.coeffs(2:2:2*W.N,:) ;
        Wr = sqrt(roots(Dw)) ;
    case 2 % near w0 approximation (beta<<1)
        BETA = math.poly.X(50) ;
        W2 = 3*xi*a^2*inv1mP(x*xi*BETA*BETA) ;
        Db = sinc_sqrt(x*W2*(alpha+2*BETA)).*(1+BETA) - m*sinc_sqrt(x*W2*(alpha-2*BETA)).*(1-BETA) ;
        Db.coeffs = Db.coeffs(1:2:2*BETA.N+1,:) ;
        betaR = roots(Db).^.5 ;
        Wr = sqrt(3*xi*a^2./(1-xi*x*(betaR).^2)) ;
    case 3 % very high frequency approximation
        n = 1:30 ;
        Wr = [.25*(4*n+m+1)*pi.*1/sqrt(x) .25*(4*n-m+1)*pi.*sqrt(xi)] ;
end

wr = (2/b)*sqrt(G/rho)*Wr ;
tol = (eps) ; wr(abs(imag(wr/2/pi/fc))>tol) = [] ;
% f = (.5*b*wt./vl).^2.*beta ; P = Kl2Beta ; 
% clf ; 
% plot(abs(beta),abs(f))
% drawnow ; set(gca,'ylimmode','manual')
% plot(abs(beta),abs(P(beta)),':')

clf ; 
plot(wt/2/pi/fc,log10(abs(d0)),'-') ;
plot(wr(:)'/2/pi/fc.*[1;1],get(gca,'ylim'),'-k','linewidth',.5) ;
% plot(((1:1000)*pi*vl/b)/2/pi/fc.*[1;1],get(gca,'ylim'),'k','linewidth',.5) ;
% plot(((1:1000)*pi*sqrt(xi*G/rho)/b)/2/pi/fc.*[1;1],get(gca,'ylim'),'b','linewidth',.5) ;
% 
% w = (n*pi*vl/b) = (n*pi*vinf/b)
% W = .5*b*w/vt = .5*n*pi/sqrt(x) = .5*n*pi*sqrt(xi)
return

%% DETERMINANTS

[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(geo,mat) ;
a = b/h ; % section aspect ratio
x = G./Q ;
alpha = 1 + 1./(x.*xi) ;
vt = sqrt(G./rho) ;

% Adimensionnal frequency & wavenumber
W = .5*h*Wd./vt ; % W2 = (.5*h*w/vt)
B = .5*h*Kd ;

% Trigonometric function approximations
sinc_sqrt = @(B)sin(sqrt(B))./sqrt(B) ;
cos_sqrt = @(B)cos(sqrt(B)) ;

% Frequency
B2 = B.*B ;
W2 = W.*W ;
% Adim. frequencies
Wl2 = x.*W2 ;
Wt2 = W2 ;
% Wavenumbers
Bl2 = Wl2 - B2 ;
Bt2 = Wt2 - B2 ;
Bst2 = W2 - B2 - 3.*xi ;
sc = 12.*xi.*(alpha-1) ; % scaling factor
sqrtAlpha2m4Beta2 = sqrt(sc).*sqrt(1+W2*((alpha-2).^2./sc)) ;
Bb2 = .5*x.*W.*(alpha.*W + sqrtAlpha2m4Beta2) - B2 ;
Bsa2 = .5*x.*W.*(alpha.*W - sqrtAlpha2m4Beta2) - B2 ;
% Trigonometric functions
cst = cos_sqrt(a^2.*Bst2) ; sst = sinc_sqrt(a^2.*Bst2) ; 
cb = cos_sqrt(a^2.*Bb2) ; sb = sinc_sqrt(a^2.*Bb2) ; 
csa = cos_sqrt(a^2.*Bsa2) ; ssa = sinc_sqrt(a^2.*Bsa2) ; 
% Determinants
Dops = 4.*x.*B2.*Bb2.*Bsa2.*(Bb2-Bsa2).*sb.*ssa.*cst + sst.*(...
            + Bsa2.*(Bb2+(1-2.*x).*B2).*(2.*x.*B2.*(Bt2-Bst2) + (Bl2-Bsa2).*(Bst2-B2)).*ssa.*cb ...
            - Bb2.*(Bsa2+(1-2.*x).*B2).*(2.*x.*B2.*(Bt2-Bst2) + (Bl2-Bb2).*(Bst2-B2)).*sb.*csa ...
        ) ; % symetric modes
Dopa = 4.*x.*B2.*Bst2.*(Bb2-Bsa2).*sst.*cb.*csa + cst.*( ...
             + (Bb2+(1-2.*x).*B2).*(2.*x.*B2.*(Bt2-Bst2) + (Bl2-Bsa2).*(Bst2-B2)).*sb.*csa ...
             - (Bsa2+(1-2.*x).*B2).*(2.*x.*B2.*(Bt2-Bst2) + (Bl2-Bb2).*(Bst2-B2)).*ssa.*cb ...
        ) ; % antisymetric modes
 

Bl = sqrt(Bl2) ; 
Bt = sqrt(Bt2) ;
Bst = sqrt(Bst2) ;
Bb = sqrt(Bb2) ; 
Bsa = sqrt(Bsa2) ;
sinc = @(x)sin(x)./x ;

D1 = 4.*x.*B.^2.*Bst.*Bb.*Bsa.*(Bb.^2-Bsa.^2).*sin(a*Bst).*cos(a*Bb).*cos(a*Bsa) ...
    + cos(a*Bst).*( ...
         + Bsa.*(Bb.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bsa.^2).*(Bst.^2-B.^2)).*sin(a*Bb).*cos(a*Bsa) ...
         - Bb.*(Bsa.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bb.^2).*(Bst.^2-B.^2)).*sin(a*Bsa).*cos(a*Bb) ...
     ) ; 
 
D2 = 4.*x.*B.^2.*Bst.^2.*Bb.*Bsa.*(Bb.^2-Bsa.^2).*sin(a*Bst).*(cos(a*Bb+a*Bsa) + cos(a*Bb-a*Bsa)) ... 2 cos(Bb).*cos(Bsa) = cos(Bb+Bsa) + cos(Bb-Bsa)
    + Bst.*cos(a*Bst).*( ...
         + Bsa.*(Bb.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bsa.^2).*(Bst.^2-B.^2)).*(sin(a*Bb+a*Bsa) + sin(a*Bb-a*Bsa)) ... 2 sin(Bb).*cos(Bsa) = sin(Bb+Bsa) + sin(Bb-Bsa)
         - Bb.*(Bsa.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bb.^2).*(Bst.^2-B.^2)).*(sin(a*Bb+a*Bsa) - sin(a*Bb-a*Bsa)) ... 2 sin(Bsa).*cos(Bb) = sin(Bb+Bsa) - sin(Bb-Bsa)
     ); 
 
D4 = 4.*x.*B.^2.*Bst.^2.*Bb.*Bsa.*(Bb.^2-Bsa.^2).*sin(a*Bst).*(cos(a*Bb+a*Bsa) + cos(a*Bb-a*Bsa)) ... 
    + Bst.*cos(a*Bst).*( ...
         + (sin(a*Bb+a*Bsa).*(Bb-Bsa) - sin(a*Bb-a*Bsa).*(Bb+Bsa)).*(...
            + ( Bb.^2.*Bsa.^2 + (Bb.^2+Bsa.^2).*(1-2.*x).*B.^2 ).*(Bst.^2-B.^2) ...
            + B.^2*(1-2.*x).*( (B.^2-Bst.^2).*Bl.^2 + 2.*x.*B.^2.*(Bst.^2-Bt.^2) ) ...
         ) ...
         + (sin(a*Bb+a*Bsa).*(Bb-Bsa) + sin(a*Bb-a*Bsa).*(Bb+Bsa)).*(...
                Bb.*Bsa.*( (Bl.^2+(1-2.*x).*B.^2).*(Bst.^2-B.^2) + 2.*x.*B.^2.*(Bt.^2-Bst.^2) ) ...
            ) ...
     ) ;
 
D5 = 4.*x.*B.^2.*Bst.^2.*Bb.*Bsa.*sinc(a*Bst).*(cos(a*Bb+a*Bsa)+cos(a*Bb-a*Bsa)) ... 
    + cos(a*Bst).*( ...
         + (sinc(a*Bb+a*Bsa) - sinc(a*Bb-a*Bsa)).*(...
            + ( Bb.^2.*Bsa.^2 + (Bb.^2+Bsa.^2).*(1-2.*x).*B.^2 ).*(Bst.^2+B.^2) ...
            + B.^2.*(1-2.*x).*(2.*x.*B.^2 - Bl.^2).*(Bst.^2+B.^2) ...
            - 2.*B.^2.*( Bb.^2.*Bsa.^2 + (Bb.^2+Bsa.^2).*(1-2.*x).*B.^2 ) ...
            + B.^2.*(1-2.*x).*(- 2.*x.*B.^2.*(B.^2+Bt.^2) + 2*B.^2.*Bl.^2) ...
         ) ...
         + (sinc(a*Bb+a*Bsa) + sinc(a*Bb-a*Bsa)).*(...
                Bb.*Bsa.*( (Bl.^2+(1-2.*x).*B.^2).*(Bst.^2+B.^2-2*B.^2) - 2.*x.*B.^2.*(Bst.^2-Bt.^2) ) ...
            ) ...
     ) ;
 
d1 = ...
    + ( Bb.^2.*Bsa.^2 + (Bb.^2+Bsa.^2).*(1-2.*x).*B.^2 ).*(Bst.^2-B.^2) ...
    + B.^2.*(1-2.*x).*( (B.^2-Bst.^2).*Bl.^2 + 2.*x.*B.^2.*(Bst.^2-Bt.^2) ) ...
    ;
d2 = 4.*x.*B.^2.*Bst.*Bb.*Bsa.*(Bb.^2-Bsa.^2).*sin(a*Bst).*cos(a*Bb).*cos(a*Bsa) ;
d3 = cos(a*Bst).*( ...
         + Bsa.*(Bb.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bsa.^2).*(Bst.^2-B.^2)).*sin(a*Bb).*cos(a*Bsa) ...
         - Bb.*(Bsa.^2+(1-2.*x).*B.^2).*(2.*x.*B.^2.*(Bt.^2-Bst.^2) + (Bl.^2-Bb.^2).*(Bst.^2-B.^2)).*sin(a*Bsa).*cos(a*Bb) ...
     ) ; 

D = Dopa;%./Bb./Bsa./(Bst2+B2) ; %
clf ; axis tight
surf((real(B)+imag(B)),W/2/pi*4,log10(abs(D))) ; 
shading interp
%% TEST LONGWAVES


m = ribbon.model(geo,mat,Kd,Wd) ;
V = Wd./Kd./sqrt(m.vt2) ;
x = m.vt2./m.vl2 ;
xi = analytical.rayleigh(mat.nu)^2 ;
s = sqrt(12*xi)/geo.h ;
z = sqrt(12*x)/geo.h ;
g = .5/xi ;
a = geo.b/geo.h ;
nu = mat.nu ;
w2 = geo.b^2/12 ;

k = Kd ;
kt2 = k.^2.*(V.^2-1) ; m.kt.^2 ;
kl2 = k.^2.*(x.*V.^2-1) ; m.kl.^2 ; 
kst2 = k.^2.*(V.^2-1) - s^2 ; m.kst.^2 ; 
kb2 = z*k.*V - k.^2 + g.*k.^2.*V.^2 ; m.kb.^2 ;
ksa2 = -z*k.*V - k.^2 + g.*k.^2.*V.^2 ; m.ksa.^2 ;

tancK = @(k2) ... tan(.5*geo.b*sqrt(k2))./sqrt(k2) ...
              ... (1 + (.5*geo.b)^2.*k2/3) ...
              ... 1 + w2.*k2 ...
              (1 + w2*k2 + (6/5)*(w2*k2).^2)...(2/15)*((.5*geo.b).^2.*k2).^2) ...
              ; %  approximation of tan(.5*b*k)/k
tanK = @(k2) tan(.5*geo.b*sqrt(k2)) ; %1 + w2.*k2 + 0*(6/5)*w2.^2*k2.^2 ; % aoproximation of tan(.5*b*k)/k
sincK = @(k2) ... sin(.5*geo.b*sqrt(k2))./sqrt(k2) ...
              (1 - w2*k2/2 + (3/40)*(w2*k2).^2) ...
              ;
cosK = @(k2) ... cos(.5*geo.b*sqrt(k2)) ...
              (1 - (3/2)*w2*k2 + (3/8)*(w2*k2).^2) ...
            ;

% Dops = 4.*x.*k.^2.*kb2.*ksa2.*(kb2-ksa2).*tancK(kb2).*tancK(ksa2) ...
%         + ksa2.*(kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*tancK(kst2).*tancK(ksa2) ...
%         - kb2.*(ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*tancK(kst2).*tancK(kb2) ...
%         ; % symetric modes
% 
% Dopa = 4.*x.*k.^2.*kst2.*(kb2-ksa2).*tancK(kst2) ...
%         + (kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*tancK(kb2) ...
%         - (ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*tancK(ksa2) ...
%         ; % antisymetric modes

Dops = 4.*x.*k.^2.*kb2.*ksa2.*(kb2-ksa2).*sincK(kb2).*sincK(ksa2).*cosK(kst2) ...
        + ksa2.*(kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*sincK(kst2).*sincK(ksa2).*cosK(kb2) ...
        - kb2.*(ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*sincK(kst2).*sincK(kb2).*cosK(ksa2) ...
        ; % symetric modes
%     
Dopa = 4.*x.*k.^2.*kst2.*(kb2-ksa2).*sincK(kst2).*cosK(kb2).*cosK(ksa2) ...
        + (kb2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-ksa2).*(kst2-k.^2)).*sincK(kb2).*cosK(ksa2).*cosK(kst2) ...
        - (ksa2+(1-2.*x).*k.^2).*(2.*x.*k.^2.*(kt2-kst2) + (kl2-kb2).*(kst2-k.^2)).*sincK(ksa2).*cosK(kb2).*cosK(kst2) ...
        ; % antisymetric modes

% B = w2.*Kd.^2 ;
% Dfun = @(B,V)(...
%         B*(2560*V^4*x + 960*V^4*a^4 + 10240*a^2*x*xi + 3840*V^4*a^2*x + 768*V^4*a^4*x + 240*V^4*a^6*xi - 10240*a^2*x^2*xi + 24576*a^4*x*xi^2 + 5376*a^6*x*xi^3 + (640*V^4*a^2)/xi - 15360*a^4*x^2*xi^2 - 3840*a^6*x^2*xi^3 + 5120*V^2*a^2*x^2*xi - 18816*V^2*a^4*x*xi^2 + 12288*V^2*a^4*x^2*xi - 8832*V^2*a^6*x*xi^2 + 2880*V^4*a^4*x*xi^2 + 768*V^4*a^4*x^2*xi - 1920*V^2*a^6*x*xi^3 + 3360*V^4*a^6*x*xi^2 + 192*V^4*a^6*x^2*xi - 576*V^2*a^8*x*xi^3 + 288*V^4*a^8*x*xi^2 + 7680*V^2*a^4*x^2*xi^2 + 3072*V^2*a^6*x^2*xi^2 + 1920*V^2*a^6*x^2*xi^3 + 1152*V^4*a^6*x^2*xi^2 + 288*V^4*a^8*x^2*xi^2 + 288*V^4*a^8*x^2*xi^3 + 72*V^4*a^10*x^2*xi^3 - 29440*V^2*a^2*x*xi - 20736*V^2*a^4*x*xi + 7680*V^4*a^2*x*xi + 9600*V^4*a^4*x*xi + 1152*V^4*a^6*x*xi) - B^2*(10240*V^2*x - 19200*V^4*x + 3840*V^6*x - 1536*V^4*a^4 + 720*V^6*a^4 + 72*V^6*a^6 - 10240*V^2*x^2 + 5120*V^4*x^2 - (1280*V^4)/xi^2 + (640*V^6)/xi^2 - 21504*a^2*x*xi + 15360*V^2*a^2*x - 32640*V^4*a^2*x - 4416*V^4*a^4*x + 7680*V^6*a^2*x + 2880*V^6*a^4*x + 144*V^6*a^6*x - 144*V^4*a^6*xi + (5120*V^4*x)/xi + (1280*V^6*x)/xi + 21504*a^2*x^2*xi - 19200*a^4*x*xi^2 - 2304*a^6*x*xi^3 - 15360*V^2*a^2*x^2 + 13824*V^4*a^2*x^2 + 3840*V^4*a^4*x^2 + 768*V^6*a^2*x^2 + 480*V^6*a^4*x^2 - (3264*V^4*a^2)/xi + (1920*V^6*a^2)/xi + (192*V^6*a^2)/xi^2 + (288*V^6*a^4)/xi + 16896*a^4*x^2*xi^2 + 2304*a^6*x^2*xi^3 - 41472*V^2*a^2*x^2*xi + (3072*V^4*a^2*x)/xi + 24576*V^2*a^4*x*xi^2 - 19200*V^2*a^4*x^2*xi + 15360*V^4*a^2*x^2*xi - 768*V^2*a^4*x^3*xi + (1152*V^6*a^2*x)/xi + 4896*V^2*a^6*x*xi^2 - 5760*V^4*a^4*x*xi^2 + 6528*V^4*a^4*x^2*xi + 1152*V^2*a^6*x*xi^3 + 384*V^4*a^4*x^3*xi + (96*V^6*a^4*x)/xi - 2592*V^4*a^6*x*xi^2 - 288*V^4*a^6*x^2*xi + 2304*V^6*a^4*x^2*xi + 72*V^2*a^8*x*xi^3 + 1152*V^4*a^6*x^3*xi - 72*V^4*a^8*x*xi^2 + 1008*V^6*a^6*x^2*xi + 72*V^6*a^6*x^3*xi - 19968*V^2*a^4*x^2*xi^2 - 2304*V^2*a^6*x^2*xi^2 + 5760*V^4*a^4*x^2*xi^2 - 1152*V^2*a^6*x^2*xi^3 - 1152*V^2*a^6*x^3*xi^2 - 1152*V^4*a^6*x^2*xi^2 + 288*V^2*a^8*x^2*xi^3 + 576*V^4*a^6*x^3*xi^2 - 288*V^2*a^8*x^3*xi^3 - 432*V^4*a^8*x^2*xi^2 + 864*V^6*a^6*x^2*xi^2 - 288*V^4*a^8*x^2*xi^3 + 288*V^4*a^8*x^3*xi^2 + 144*V^4*a^8*x^3*xi^3 + 324*V^6*a^8*x^2*xi^2 + 108*V^6*a^8*x^3*xi^2 + 27*V^6*a^10*x^3*xi^3 + 60864*V^2*a^2*x*xi + 24960*V^2*a^4*x*xi - 28416*V^4*a^2*x*xi - 24768*V^4*a^4*x*xi + 2880*V^6*a^2*x*xi - 2016*V^4*a^6*x*xi + 4320*V^6*a^4*x*xi + 1008*V^6*a^6*x*xi + 36*V^6*a^8*x*xi) + 10240*a^4*x*xi^2 + 3072*a^6*x*xi^3 - 3840*V^2*a^4*x*xi^2 - 3840*V^2*a^6*x*xi^2 - 960*V^2*a^6*x*xi^3 - 960*V^2*a^8*x*xi^3 - 2560*V^2*a^2*x*xi - 2560*V^2*a^4*x*xi ...
%     ) ;
% Dapp = arrayfun(Dfun,B+V*0,V+B*0) ;

% fm = sqrt(m.vt2) + 0*Wd ; f = sqrt(mat.G/mat.rho) + 0*Wd ;
% fm = m.kt.^2 ; f = kt2 ;
% fm = m.kl.^2 ; f = kl2 ;
% fm = m.kst.^2 ; f = kst2 ;
% fm = m.kb.^2 ; f = kb2 ;
% fm = m.ksa.^2 ; f = ksa2 ;
%k2 = ksa2 ; fm = tancK(k2) ; f = (1 + ((.5*geo.b).^2.*k2)/3 + (2/15)*((.5*geo.b).^2.*k2).^2)*(.5*geo.b) ;
f = Dopa ; fm = m.Dopa*2*x./m.kb./m.ksa ;
%fm = Dopa ; f = Dapp ;
%f = Dops ; fm = m.Dops*2*x ;

plotFun = @(f)log10(abs(f)) ;
rootFun = @(D)abs(D)./medfilt2(abs(D),5.*[1 1],'symmetric') ;
Dmin = .3 ; Dmax = 1 ; 
scaleFun = @(D)(max(Dmin,min(D,Dmax))-Dmin)./(Dmax-Dmin) ;

clf ;
%%%%%
mysubplot(1,3,1) ; axis tight ; title('model')
surf(real(Kd)+imag(Kd),real(Wd),plotFun(fm),'facecolor','interp','edgecolor','none')
mysubplot(1,3,2) ; axis tight ; title('approximation')
surf(real(Kd)+imag(Kd),real(Wd),plotFun(f),'facecolor','interp','edgecolor','none')
mysubplot(1,3,3) ; axis tight ; title('roots')
surf(real(Kd)+imag(Kd),real(Wd),0*abs(f),1 - (1-scaleFun(rootFun(fm))).*reshape(1-[1 0 0],1,1,3) - (1-scaleFun(rootFun(f))).*reshape(1-[0 0 1],1,1,3),'facecolor','interp','edgecolor','none')
linkprop(findobj(gcf,'type','axes'),'view') ;
%%%%%
% axis tight ;
% surf(real(Kd)+imag(Kd),real(Wd),plotFun(fm),'facecolor','interp','edgecolor','none')
% surf(real(Kd)+imag(Kd),real(Wd),plotFun(f),'facecolor','none','edgecolor','k')
%%%%%
%view([50 30])
B = w2*ks.^2 ;
Vr2 = ...8*xi*a.^2./(a.^2+1)./(3*xi.*a.^2+2) ...
      ... 2*ks.^2.*geo.h^2*(nu + 1)/12 ...
      (16*a^2*xi*(3*xi*a^2 + 10))/(5*(a^2 + 1)*(3*a^4*xi^2 + 12*a^2*xi + 8)) ...
      ...+ 1*B.*(-(4*(64000*x + 128000*a^2*x + 64000*a^4*x + 134400*a^2*xi + 345600*a^4*xi + 211200*a^6*xi + 64000*a^4 + 187200*a^4*xi^2 + 347520*a^6*xi^2 + 96480*a^6*xi^3 + 160320*a^8*xi^2 + 37440*a^8*xi^3 + 32760*a^8*xi^4 - 59040*a^10*xi^3 - 60480*a^10*xi^4 + 7200*a^10*xi^5 - 93240*a^12*xi^4 - 23040*a^12*xi^5 + 675*a^12*xi^6 - 30240*a^14*xi^5 - 2430*a^14*xi^6 - 3105*a^16*xi^6 + 160000*a^2*x*xi + 140800*a^4*x*xi - 19200*a^6*x*xi + 81600*a^4*x*xi^2 - 120960*a^6*x*xi^2 - 67200*a^6*x*xi^3 - 144960*a^8*x*xi^2 - 188160*a^8*x*xi^3 - 70200*a^8*x*xi^4 - 83952*a^10*x*xi^4 - 20700*a^10*x*xi^5 + 64872*a^12*x*xi^4 - 16488*a^12*x*xi^5 - 2025*a^12*x*xi^6 + 24948*a^14*x*xi^5 - 1242*a^14*x*xi^6 + 2727*a^16*x*xi^6 - 64000))/(125*(a^2 + 1)^3*(3*a^4*xi^2 + 12*a^2*xi + 8)^3)) ...
       + 1*B.*((64*(9*a^10*x^2*xi^4 + 36*a^8*x^2*xi^4 + 36*a^8*x^2*xi^3 + 36*a^8*x*xi^3 + 144*a^6*x^2*xi^3 + 24*a^6*x^2*xi^2 + 420*a^6*x*xi^3 + 144*a^6*x*xi^2 + 30*a^6*xi^2 + 96*a^4*x^2*xi^2 + 360*a^4*x*xi^3 + 1200*a^4*x*xi^2 + 96*a^4*x*xi + 120*a^4*xi + 960*a^2*x*xi^2 + 480*a^2*x*xi + 80*a^2 + 320*x*xi)*(- 2727*a^16*x^2*xi^6 + 3105*a^16*x*xi^6 + 1242*a^14*x^2*xi^6 - 24948*a^14*x^2*xi^5 + 2430*a^14*x*xi^6 + 25272*a^14*x*xi^5 + 2025*a^12*x^2*xi^6 + 16488*a^12*x^2*xi^5 - 64872*a^12*x^2*xi^4 - 675*a^12*x*xi^6 + 19800*a^12*x*xi^5 + 43848*a^12*x*xi^4 + 2160*a^12*xi^4 + 20700*a^10*x^2*xi^5 + 83952*a^10*x^2*xi^4 - 7200*a^10*x*xi^5 + 29520*a^10*x*xi^4 - 101088*a^10*x*xi^3 + 23040*a^10*xi^3 + 70200*a^8*x^2*xi^4 + 188160*a^8*x^2*xi^3 + 144960*a^8*x^2*xi^2 - 32760*a^8*x*xi^4 - 127680*a^8*x*xi^3 - 306240*a^8*x*xi^2 + 87360*a^8*xi^2 + 67200*a^6*x^2*xi^3 + 120960*a^6*x^2*xi^2 + 19200*a^6*x^2*xi - 96480*a^6*x*xi^3 - 385920*a^6*x*xi^2 - 96000*a^6*x*xi + 134400*a^6*xi - 81600*a^4*x^2*xi^2 - 140800*a^4*x^2*xi - 64000*a^4*x^2 - 187200*a^4*x*xi^2 - 179200*a^4*x*xi + 64000*a^4*x + 64000*a^4 - 160000*a^2*x^2*xi - 128000*a^2*x^2 - 134400*a^2*x*xi + 128000*a^2*x - 64000*x^2 + 64000*x))/(125*x*xi*(a^2 + 1)^3*(3*a^4*xi^2 + 12*a^2*xi + 8)^3*(5120*x + 7680*a^2*x + 1536*a^4*x + 480*a^6*xi + 1920*a^4 + (1280*a^2)/xi + 15360*a^2*x*xi + 19200*a^4*x*xi + 2304*a^6*x*xi + 5760*a^4*x*xi^2 + 1536*a^4*x^2*xi + 6720*a^6*x*xi^2 + 384*a^6*x^2*xi + 576*a^8*x*xi^2 + 2304*a^6*x^2*xi^2 + 576*a^8*x^2*xi^2 + 576*a^8*x^2*xi^3 + 144*a^10*x^2*xi^3))) ...
      ;
vr = sqrt(Vr2*(mat.G/mat.rho)) ;
wr = ks.*vr ;
plot(real(ks)+imag(ks),real(wr),':k')
%%
clf ; 
plot(real(ks)+imag(ks),abs(tan(.5*geo.b*ks)./ks))
plot(real(ks)+imag(ks),abs(1+(.5*geo.b)^2*ks.^2/3 + (.5*geo.b)^4*ks.^4*2/15)*(.5*geo.b))
set(gca,'yscale','log')

%% TEST TWIST MOTION COEFFICIENTS
a = linspace(1,20,1000)' ;
nu = linspace(-1,.5,10) ;
xi = 10/10 ; pi.^2./12
x = (1-nu)./2 ;
phi = (16.*a.^2.*xi.*(3.*xi.*a.^2 + 10))./(5.*(a.^2 + 1).*(3.*a.^4.*xi.^2 + 12.*a.^2.*xi + 8)) ;
psi = ...(-(4.*(64000.*x + 128000.*a.^2.*x + 64000.*a.^4.*x + 134400.*a.^2.*xi + 345600.*a.^4.*xi + 211200.*a.^6.*xi + 64000.*a.^4 + 187200.*a.^4.*xi.^2 + 347520.*a.^6.*xi.^2 + 96480.*a.^6.*xi.^3 + 160320.*a.^8.*xi.^2 + 37440.*a.^8.*xi.^3 + 32760.*a.^8.*xi.^4 - 59040.*a.^10.*xi.^3 - 60480.*a.^10.*xi.^4 + 7200.*a.^10.*xi.^5 - 93240.*a.^12.*xi.^4 - 23040.*a.^12.*xi.^5 + 675.*a.^12.*xi.^6 - 30240.*a.^14.*xi.^5 - 2430.*a.^14.*xi.^6 - 3105.*a.^16.*xi.^6 + 160000.*a.^2.*x.*xi + 140800.*a.^4.*x.*xi - 19200.*a.^6.*x.*xi + 81600.*a.^4.*x.*xi.^2 - 120960.*a.^6.*x.*xi.^2 - 67200.*a.^6.*x.*xi.^3 - 144960.*a.^8.*x.*xi.^2 - 188160.*a.^8.*x.*xi.^3 - 70200.*a.^8.*x.*xi.^4 - 83952.*a.^10.*x.*xi.^4 - 20700.*a.^10.*x.*xi.^5 + 64872.*a.^12.*x.*xi.^4 - 16488.*a.^12.*x.*xi.^5 - 2025.*a.^12.*x.*xi.^6 + 24948.*a.^14.*x.*xi.^5 - 1242.*a.^14.*x.*xi.^6 + 2727.*a.^16.*x.*xi.^6 - 64000))./(125.*(a.^2 + 1).^3.*(3.*a.^4.*xi.^2 + 12.*a.^2.*xi + 8).^3)) ;
      ...(64.*(9.*a.^10.*x.^2.*xi.^4 + 36.*a.^8.*x.^2.*xi.^4 + 36.*a.^8.*x.^2.*xi.^3 + 36.*a.^8.*x.*xi.^3 + 144.*a.^6.*x.^2.*xi.^3 + 24.*a.^6.*x.^2.*xi.^2 + 420.*a.^6.*x.*xi.^3 + 144.*a.^6.*x.*xi.^2 + 30.*a.^6.*xi.^2 + 96.*a.^4.*x.^2.*xi.^2 + 360.*a.^4.*x.*xi.^3 + 1200.*a.^4.*x.*xi.^2 + 96.*a.^4.*x.*xi + 120.*a.^4.*xi + 960.*a.^2.*x.*xi.^2 + 480.*a.^2.*x.*xi + 80.*a.^2 + 320.*x.*xi).*(- 2727.*a.^16.*x.^2.*xi.^6 + 3105.*a.^16.*x.*xi.^6 + 1242.*a.^14.*x.^2.*xi.^6 - 24948.*a.^14.*x.^2.*xi.^5 + 2430.*a.^14.*x.*xi.^6 + 25272.*a.^14.*x.*xi.^5 + 2025.*a.^12.*x.^2.*xi.^6 + 16488.*a.^12.*x.^2.*xi.^5 - 64872.*a.^12.*x.^2.*xi.^4 - 675.*a.^12.*x.*xi.^6 + 19800.*a.^12.*x.*xi.^5 + 43848.*a.^12.*x.*xi.^4 + 2160.*a.^12.*xi.^4 + 20700.*a.^10.*x.^2.*xi.^5 + 83952.*a.^10.*x.^2.*xi.^4 - 7200.*a.^10.*x.*xi.^5 + 29520.*a.^10.*x.*xi.^4 - 101088.*a.^10.*x.*xi.^3 + 23040.*a.^10.*xi.^3 + 70200.*a.^8.*x.^2.*xi.^4 + 188160.*a.^8.*x.^2.*xi.^3 + 144960.*a.^8.*x.^2.*xi.^2 - 32760.*a.^8.*x.*xi.^4 - 127680.*a.^8.*x.*xi.^3 - 306240.*a.^8.*x.*xi.^2 + 87360.*a.^8.*xi.^2 + 67200.*a.^6.*x.^2.*xi.^3 + 120960.*a.^6.*x.^2.*xi.^2 + 19200.*a.^6.*x.^2.*xi - 96480.*a.^6.*x.*xi.^3 - 385920.*a.^6.*x.*xi.^2 - 96000.*a.^6.*x.*xi + 134400.*a.^6.*xi - 81600.*a.^4.*x.^2.*xi.^2 - 140800.*a.^4.*x.^2.*xi - 64000.*a.^4.*x.^2 - 187200.*a.^4.*x.*xi.^2 - 179200.*a.^4.*x.*xi + 64000.*a.^4.*x + 64000.*a.^4 - 160000.*a.^2.*x.^2.*xi - 128000.*a.^2.*x.^2 - 134400.*a.^2.*x.*xi + 128000.*a.^2.*x - 64000.*x.^2 + 64000.*x))./(125.*x.*xi.*(a.^2 + 1).^3.*(3.*a.^4.*xi.^2 + 12.*a.^2.*xi + 8).^3.*(5120.*x + 7680.*a.^2.*x + 1536.*a.^4.*x + 480.*a.^6.*xi + 1920.*a.^4 + (1280.*a.^2)./xi + 15360.*a.^2.*x.*xi + 19200.*a.^4.*x.*xi + 2304.*a.^6.*x.*xi + 5760.*a.^4.*x.*xi.^2 + 1536.*a.^4.*x.^2.*xi + 6720.*a.^6.*x.*xi.^2 + 384.*a.^6.*x.^2.*xi + 576.*a.^8.*x.*xi.^2 + 2304.*a.^6.*x.^2.*xi.^2 + 576.*a.^8.*x.^2.*xi.^2 + 576.*a.^8.*x.^2.*xi.^3 + 144.*a.^10.*x.^2.*xi.^3)) ;
    (4.*(- 2727.*a.^16.*x.^2.*xi.^6 + 3105.*a.^16.*x.*xi.^6 + 1242.*a.^14.*x.^2.*xi.^6 - 24948.*a.^14.*x.^2.*xi.^5 + 2430.*a.^14.*x.*xi.^6 + 25272.*a.^14.*x.*xi.^5 + 2025.*a.^12.*x.^2.*xi.^6 + 16488.*a.^12.*x.^2.*xi.^5 - 64872.*a.^12.*x.^2.*xi.^4 - 675.*a.^12.*x.*xi.^6 + 19800.*a.^12.*x.*xi.^5 + 43848.*a.^12.*x.*xi.^4 + 2160.*a.^12.*xi.^4 + 20700.*a.^10.*x.^2.*xi.^5 + 83952.*a.^10.*x.^2.*xi.^4 - 7200.*a.^10.*x.*xi.^5 + 29520.*a.^10.*x.*xi.^4 - 101088.*a.^10.*x.*xi.^3 + 23040.*a.^10.*xi.^3 + 70200.*a.^8.*x.^2.*xi.^4 + 188160.*a.^8.*x.^2.*xi.^3 + 144960.*a.^8.*x.^2.*xi.^2 - 32760.*a.^8.*x.*xi.^4 - 127680.*a.^8.*x.*xi.^3 - 306240.*a.^8.*x.*xi.^2 + 87360.*a.^8.*xi.^2 + 67200.*a.^6.*x.^2.*xi.^3 + 120960.*a.^6.*x.^2.*xi.^2 + 19200.*a.^6.*x.^2.*xi - 96480.*a.^6.*x.*xi.^3 - 385920.*a.^6.*x.*xi.^2 - 96000.*a.^6.*x.*xi + 134400.*a.^6.*xi - 81600.*a.^4.*x.^2.*xi.^2 - 140800.*a.^4.*x.^2.*xi - 64000.*a.^4.*x.^2 - 187200.*a.^4.*x.*xi.^2 - 179200.*a.^4.*x.*xi + 64000.*a.^4.*x + 64000.*a.^4 - 160000.*a.^2.*x.^2.*xi - 128000.*a.^2.*x.^2 - 134400.*a.^2.*x.*xi + 128000.*a.^2.*x - 64000.*x.^2 + 64000.*x))./(125.*x.*(a.^2 + 1).^3.*(3.*a.^4.*xi.^2 + 12.*a.^2.*xi + 8).^3) ; 
clf ; 
plot3(a+x*0,psi,x+a*0,':')
plot3(a+x*0,phi,x+a*0,'-k')
%plot3(a+x*0,psi./phi,nu+a*0,'-.')

A = {a.^0.*nu.^0 a.^0.*nu.^1 a.^-1*nu.^0 a.^-1*nu.^1} ;
A = cat(ndims(A{1})+1,A{:}) ;
Ap = (phi+0*psi).*A ;
Ap = reshape(Ap,[],size(Ap,ndims(Ap))) ;
cc = Ap \ psi(:)
set(gca,'colororderindex',1)
plot3(a+x*0,reshape(Ap*cc,size(A(:,:,1))),x+a*0,'-','linewidth',.1)
%% WAVENUMBER COMPUTATION
    ff = linspace(0,1.1,100)*fc ;
    [Uf,Kf] = safe.computeK(mesh,mat,2*pi*ff,nModes) ;
    Ff = repmat(ff,[nModes 1]) ;
    clf ; plot3(real(Kf(:))/kc,Ff(:)/fc,imag(Kf(:))/kc,'.','displayname','SAFE')
    
%% RIBBON MODEL TRACKING
% Ribbon Model determinant
    nPix = 1000 ;
    [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,kp,2*pi*f,nPix) ;
% Display
    clf ; hold on
    surf(real(Kd)+imag(Kd),Wd,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP-S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP-A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP-S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP-A') ;
% Zero-K frequencies
    n = 0:100 ;
    wl0 = n*pi*sqrt(mat.Q./mat.rho)./geo.b ;
    wt0 = n*pi*sqrt(mat.G./mat.rho)./geo.b ;
    xi = (pi^2/12) ; % shear correction factor
    eta = 12./(geo.h.^2) ; % shape factor
    wst0 = sqrt((mat.G./mat.rho).*(4*pi^2*n.^2./geo.b.^2 + eta.*xi)) ;
    plot(n*0,real(wl0),'x')
    plot(n*0,real(wt0),'x')
    plot(n*0,real(wst0),'x')
% Legends
    xlabel 'Wavenumber $k$ (rad/mm)' ; 
    ylabel 'Frequency $\omega$ (rad/s)'
    lgd = legend ; lgd.Location = 'southeast' ;
% Rescaling
    obj = allchild(gca) ;
    set(obj,{'xdata'},cellfun(@(x)x./kc,get(obj,'XData'),'uni',false))
    set(obj,{'ydata'},cellfun(@(x)x./2/pi/fc,get(obj,'YData'),'uni',false))
    xlabel 'Normalized Wavenumber $\frac{k \times h}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2h}{c_t}$'
    set(gca,'xlim',[-.5 1],'ylim',[0 1.1]) ;
    
    
%% ASYMPTOTIC VELOCITY
clc,clear all
nPix = 1000 ;
% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*1 ; % width (mm)
    dx = min([geo.b geo.h])/10 ; % element size
% Material
    mat = [] ;
    mat.E = 70e3*(1+0.00000i) ; % Young modulus (MPa)
    mat.nu = 40/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
    %mat.xi = analytical.rayleigh(mat.nu)^2 ;
% Adimensionalized numbers
    fc = .5*sqrt(real(mat.G)/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Phase speed domain
    ct = sqrt(mat.G./mat.rho) ;
    cR = analytical.rayleigh(mat.nu) ;
    %cR = analytical.rayleigh(mat.nu/(1+mat.nu)) ;
    %cR = analytical.rayleigh(mat.nu/(1+2*mat.nu)) ;
    cE = analytical.edge(mat.nu) ;
    c = ct*linspace(.9,1.1,nPix)' ; 
% Frequency domain
    f = fc*linspace(0,50,nPix) ; 
% Ribbon Model determinant
    [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,2*pi*f./c,2*pi*f) ;
% Display
    clf ; hold on
    surf(Wd./2./pi./fc,Wd./Kd./ct,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP/S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP/A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP/S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP/A') ;
    plot(get(gca,'xlim'),[1 1]*cR,':','displayname','Vr') ;
    plot(get(gca,'xlim'),[1 1]*cE,':','displayname','Ve') ;
% Legends
    xlabel 'Normalized Frequency $\frac{f \times 2h}{c_t}$' ; 
    ylabel 'Normalized Velocity $\frac{c}{c_t}$'
    lgd = legend ; %lgd.Location = 'southeast' ;





