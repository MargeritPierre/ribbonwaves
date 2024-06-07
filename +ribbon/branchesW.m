function [k,w] = branchesW(GEO,MAT,k,n_or_wrange)
%% Dispersion diagram computation
GEO = geo ; MAT = mat ; k = ks(:).' ; n_or_wrange = uint8(4) ;

% Adimensionalized numbers
    wc = pi*sqrt(real(MAT.G)/MAT.rho)/GEO.h ; % height-cutoff frequency
    kc = pi/GEO.h ; % height-cutoff wavenumber
    kshift = 5e-4*kc ;

% Start with cutoff frequencies
    [wips0,wipa0,wops0,wopa0] = ribbon.cutoffW(GEO,MAT,n_or_wrange) ;
    [~,ik0] = min(abs(k)) ;
    if abs(k(ik0))>kshift ; kshift = k(ik0) ; 
    else ; kshift = (k(ik0)+eps)/(abs(k(ik0))+eps)*kshift ; 
    end

% Init display
    clf ; hold on
    surf((real(Kd)+imag(Kd))/kc,Wd/wc,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP-S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP-A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP-S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP-A') ;
    plot3(real(Ks(:))/kc+imag(Ks(:))/kc,real(Ws(:))/wc,imag(Ws(:))/wc,'.k','displayname','SAFE') ; 
    %legend

% In-Plane Symmetric 
kips = k + 0*wips0 ;
wips = NaN*kips ; wips(:,ik0) = wips0 ;
% The zero-frequency mode is the longitudinal mode with phase velocity vl
kips(wips0==0,ik0) = kshift ;
wips(wips0==0,ik0) = sqrt(MAT.Q./MAT.rho).*kshift ;
% Define the cost function
Fips = @(kn,wn)ribbon.model(GEO,MAT,kc*kn,wc*wn).Dips_n ; 
% Build the branches
for ww = 1:numel(wips0)
    wips(ww,:) = wc*zerocurve(Fips,kips(ww,:)/kc,wips(ww,ik0)/wc) ;
end
% Display
plot3(real(kips)'/kc+imag(kips)'/kc,real(wips)'/wc,imag(wips)'/wc,':','color',colors(1,:)) ;

% In-Plane Anti-symmetric 
kipa = k + 0*wipa0 ;
wipa = NaN*kipa ; wipa(:,ik0) = wipa0 ;
% The zero-frequency mode is the bending mode along width
kipa(wipa0==0,ik0) = kshift ;
wipa(wipa0==0,ik0) = sqrt(MAT.E*GEO.b^2./MAT.rho/12)*(kshift.^2) ;
% Define the cost function
Fipa = @(kn,wn)ribbon.model(GEO,MAT,kc*kn,wc*wn).Dipa_n ; 
% Build the branches
for ww = 1:numel(wipa0)
    wipa(ww,:) = wc*zerocurve(Fipa,kipa(ww,:)/kc,wipa(ww,ik0)/wc) ;
end
% Display
plot3(real(kipa)'/kc+imag(kipa)'/kc,real(wipa)'/wc,imag(wipa)'/wc,':','color',colors(2,:)) ;

% Out-of-Plane Symmetric 
%wops0 = wops0(1)
kops = k + 0*wops0 ;
wops = NaN*kops ; wops(:,ik0) = wops0 ;
% The zero-frequency mode is the bending mode in the height direction
kops(wops0==0,ik0) = kshift ;
wops(wops0==0,ik0) = sqrt(MAT.E*GEO.h^2./MAT.rho/12)*(kops(wops0==0,ik0).^2) ;
% Define the cost function
Fops = @(kn,wn)ribbon.model(GEO,MAT,kc*kn,wc*wn).Dops_n ; 
% Build the branches
for ww = 1:numel(wops0)
    wops(ww,:) = wc*zerocurve(Fops,kops(ww,:)/kc,wops(ww,ik0)/wc) ;
end
% Display
plot3(real(kops)'/kc+imag(kops)'/kc,real(wops)'/wc,imag(wops)'/wc,':','color',colors(3,:)) ;

% Out-of-Plane Anti-symmetric 
kopa = k + 0*wopa0 ;
wopa = NaN*kopa ; wopa(:,ik0) = wopa0 ;
% The zero-frequency mode is the twist mode 
kopa(wopa0==0,ik0) = kshift ;
wopa(wopa0==0,ik0) = sqrt(MAT.G./MAT.rho)*kshift ;
% Define the cost function
Fopa = @(kn,wn)ribbon.model(GEO,MAT,kc*kn,wc*wn).Dopa_n ; 
% Build the branches
for ww = 1:numel(wopa0)
    wopa(ww,:) = wc*zerocurve(Fopa,kopa(ww,:)/kc,wopa(ww,ik0)/wc) ;
end
% Display
plot3(real(kopa)'/kc+imag(kopa)'/kc,real(wopa)'/wc,imag(wopa)'/wc,':','color',colors(4,:)) ;


end

function y = zerocurve(F,x,y0)
% Build the parametrized curve y(x) so that F(x,y(x))=0
% y0 is the initial "guess" for y(x=0)
    y = NaN(size(x)) ;
% Use the guess at x==0
    [~,i0] = min(abs(x)) ;
    y(i0) = findzero(@(y)F(x(i0),y),y0) ;
% Go increasing from x0
    for ii = i0+1:numel(x)
        y(ii) = y(ii-1) ;
        y(ii) = y(ii) + dy_dx(F,x(ii-1),y(ii-1))*(x(ii)-x(ii-1)) ;
        y(ii) = findzero(@(y)F(x(ii),y),y(ii)) ;
    end
% Go decreasing from x0
    for ii = i0-1:-1:1
        y(ii) = y(ii+1) ;
        y(ii) = y(ii) + dy_dx(F,x(ii+1),y(ii+1))*(x(ii)-x(ii+1)) ;
        y(ii) = findzero(@(y)F(x(ii),y),y(ii)) ;
    end
end

function d = dy_dx(F,x0,y0)
% uses the fact that F(x,y) = 0 
%   thus dF = dF_dx.dx + dF_dy.dy = 0 
%   to compute dy/dx
    [~,dF_dy] = feval(@(y)F(x0,y),y0) ;
    [~,dF_dx] = feval(@(x)F(x,y0),x0) ;
    d = - (dF_dx/dF_dy) ;
end

function p = findzero(F,p,tol)
% Find the zero of function F(p) closest to p
if nargin<3 ; tol = 1e-6 ; end
    dp = inf ; lmbda = 1 ; r = 1.5 ;
    [f0,j0] = feval(F,p) ;
    while norm(dp)>tol
        dp = -lmbda*conj(j0)*f0/norm(j0)^2 ;
        [f,j] = feval(F,p+dp) ;
        if abs(f0)<abs(f) % force decreasing function
            %return ;
            lmbda = lmbda/r ;
        else
            lmbda = min(lmbda*r,inf) ;
            p = p + dp ;  
            f0 = f ; j0 = j ;
        end
    end
end

function [f,j] = feval(F,p,delta)
% Value and gradient of the function F(p) evaluated at p
if nargin<3 ; delta = 1e-6 ; end % finite difference step
    nP = numel(p) ;
    dp = [eye(nP) ; -eye(nP) ; zeros(1,nP)].*delta ;
    fp = F(p+dp) ;
    f = fp(end) ;
    j = (.5./delta).*(fp(1:nP)-fp(nP+1:2*nP)).' ;
end

