function [k,w] = branches(GEO,MAT,krange,wrange)
%% Dispersion diagram computation
GEO = geo ; MAT = mat ;

% Adimensionalized numbers
    fc = .5*sqrt(real(MAT.G)/MAT.rho)/GEO.h ; % height-cutoff frequency
    kc = pi/GEO.h ; % height-cutoff wavenumber
    
krange = [0 1]*kc ; 
wrange = [0 1.1]*2*pi*fc ; 

n = 10 ; 
prange = [krange(:)/kc wrange(:)/2/pi/fc] ;

% Start with cutoff frequencies
[wips,wipa,wops,wopa] = ribbon.cutoffW(GEO,MAT,wrange) ;

% Init display
clf ; hold on
surf((real(Kd)+imag(Kd))/kc,Wd/2/pi/fc,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
plot(NaN,NaN,'color',colors(1,:),'displayname','IP-S') ;
plot(NaN,NaN,'color',colors(2,:),'displayname','IP-A') ;
plot(NaN,NaN,'color',colors(3,:),'displayname','OP-S') ;
plot(NaN,NaN,'color',colors(4,:),'displayname','OP-A') ;
plot3(real(Ks(:))/kc+imag(Ks(:))/kc,real(Ws(:))/2/pi/fc,imag(Ws(:))/2/pi/fc,'.k','displayname','SAFE') ; 
legend
% In-Plane Symmetric 
kips = wips*0 ;
kips(wips==0) = 1e-6*kc ; 
wips(wips==0) = sqrt(MAT.Q./MAT.rho)*kips(1) ;
Fips = @(p)ribbon.model(GEO,MAT,p(:,1)*kc,2*pi*fc*p(:,2)).Dips_n ; 
pips = follow([kips(:)/kc wips(:)/2/pi/fc],Fips,prange) ;
% plot(real(pips(:,1)),real(pips(:,2)),'.') ;

% In-Plane Anti-symmetric
kipa = wipa*0 ;
kipa(wipa==0) = 20e-3*kc ; 
wipa(wipa==0) = sqrt(MAT.E*GEO.b^2./MAT.rho/12)*(kipa(1).^2) ;
Fipa = @(p)ribbon.model(GEO,MAT,p(:,1)*kc,2*pi*fc*p(:,2)).Dipa_n ; 
pipa = follow([kipa(:)/kc wipa(:)/2/pi/fc],Fipa,prange) ;
%plot(real(pipa(:,1)),real(pipa(:,2)),'.') ;

% Out-of-Plane Symmetric
kops = wops*0 ;
kops(wops==0) = 10e-3*kc ; 
wops(wops==0) = sqrt(MAT.E*GEO.h^2./MAT.rho/12)*(kops(1).^2) ;
Fops = @(p)ribbon.model(GEO,MAT,p(:,1)*kc,2*pi*fc*p(:,2)).Dops_n ; 
pops = follow([kops(:)/kc wops(:)/2/pi/fc],Fops,prange) ;
%plot(real(pipa(:,1)),real(pipa(:,2)),'.') ;

% Out-of-Plane Anti-symmetric
kopa = wopa*0 ;
kopa(wopa==0) = 10e-3*kc ; 
wopa(wopa==0) = sqrt(MAT.G./MAT.rho)*kopa(1) ;
Fopa = @(p)ribbon.model(GEO,MAT,p(:,1)*kc,2*pi*fc*p(:,2)).Dopa_n ; 
popa = follow([kopa(:)/kc wopa(:)/2/pi/fc],Fopa,prange) ;
%plot(real(pipa(:,1)),real(pipa(:,2)),'.') ;



end

function p = follow(p0,F,prange)
% Follow a set of branches
    p = cell(1,size(p0,1)) ;
    for pp = 1:size(p0,1)
        plK = plot(NaN,NaN,'o') ;
        p{pp} = zeros(0,size(p0,2)) ;
        newP = findzero(F,p0(pp,:)) ;
        while all(newP>=prange(1,:)) && all(newP<=prange(2,:))
            newP(1) = real(newP(1)) ;
            p{pp}(end+1,:) = newP ;
            newP = continuation(F,p{pp}(end,:)) ;
            plK.XData = real(p{pp}(:,1)) ; plK.YData = real(p{pp}(:,2)) ; plK.ZData = imag(p{pp}(:,2)) ; %drawnow
        end
    end
    p(end+1,:) = {NaN(1,size(p0,2))} ;
    p = cat(1,p{:}) ;
end

function p = continuation(F,p,dp,tol)
% Continue following a root with a given "guess" step dp
if nargin<4 ; dp = 1/100*[1 1] ; end
if nargin<5 ; tol = 1e-3*dp ; end
% Prediction with constant-step dp
    [~,j] = eval(F,p) ;
    t = j*[0 -1 ; 1 0] ;
    A = [j;t] ; b = [0;norm(j)*norm(dp)] ; 
    dp = sign(t*dp(:))*(A\b).' ;
    p = p + dp ;
% Correction (back to the branches)
    p = findzero(F,p,tol) ;
end

function p = findzero(F,p,tol)
% Find the zero of function F(p) closest to p
if nargin<3 ; tol = 1e-6 ; end
    dp = inf ; 
    [f0,j0] = eval(F,p) ;
    while norm(dp)>tol
        %plot(p(1),p(2),'.') ; drawnow ;
        dp = -conj(j0)*f0/norm(j0)^2 ;
        [f,j] = eval(F,p+dp) ;
        if abs(f0)<abs(f) ; return ; end % force decreasing function
        p = p + dp ; 
        f0 = f ; j0 = j ;
    end
end

function [f,j] = eval(F,p,delta)
% Value and gradient of the function F(p) evaluated at p
if nargin<3 ; delta = 1e-6 ; end % finite difference step
    dp = [1 0 ; -1 0 ; 0 1 ; 0 -1 ; 0 0]*delta ;
    fp = F(p+dp) ;
    f = fp(end) ;
    j = .5*[fp(1)-fp(2) fp(3)-fp(4)]/delta ;
end

