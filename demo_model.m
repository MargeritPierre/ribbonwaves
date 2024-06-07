%% WAVENUMBER SPECTRUM OF A PRISMATIC WAVEGUIDE
% Comparison between the SAFE method and the reduced ribbon model
clc
%clear all

% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*3 ; % width (mm)
    dx = min([geo.b geo.h])/5 ; % element size
% Material
    mat = [] ;
    mat.E = 70e3*(1+0.00000i) ; % Young modulus (MPa)
    mat.nu = 33/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
% Adimensionalized numbers
    fc = .5*sqrt(real(mat.G)/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Wavenumber domain
    kp = kc*[-.25i 0 1] ; % wavenumber keypoints
% Ribbon Model determinant
    nPix = 500 ;
    f = [0 1.1]*fc ; % frequency range
    [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,kp,2*pi*f,nPix) ;
% Eigenvalue options
    nModes = 25 ;
    dk = 0.01*kc ; 
% Build the mesh
    mesh = safe.mesh.gridmesh([geo.b geo.h],dx) ; % quad mesh
    mesh = safe.mesh.quad8(mesh) ; % serendipity element
    mesh = safe.mesh.quadrature(mesh,'quad2') ; % generate a quadrature rule
    clf ; axis equal ; safe.mesh.plotmesh(mesh) ; 
% Compute the SAFE frequencies & modes
    ks = interp1([0 cumsum(abs(diff(kp)))],kp,0:dk:sum(abs(diff(kp)))) ; % re-interpolate from keypoints
    [Us,Ws] = safe.computeW(mesh,mat,ks,nModes,false) ;
    Ks = repmat(ks,[nModes 1]) ;
% Display
    clf ; hold on
    surf(real(Kd)+imag(Kd),Wd,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP/S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP/A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP/S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP/A') ;
    plot3(real(Ks(:))+imag(Ks(:)),real(Ws(:)),imag(Ws(:)),'.k','displayname','SAFE') ; 
% Legends
    xlabel 'Wavenumber $k$ (rad/mm)' ; 
    ylabel 'Frequency $\omega$ (rad/s)'
    lgd = legend ; %lgd.Location = 'southeast' ;
% Cutoff
%     n = 0:5 ;
%     w0l = pi*sqrt(mat.Q/mat.rho)*((n(:)/geo.b).^2 + (n/geo.h).^2).^.5  ;
%     w0t = pi*sqrt(mat.G/mat.rho)*((n(:)/geo.b).^2 + (n/geo.h).^2).^.5  ;
%     plot(w0l(:)*0,w0l(:),'+','displayname','$\omega_0^\ell$') ;
%     plot(w0t(:)*0,w0t(:),'+','displayname','$\omega_0^t$') ;
% Rayleigh velocity
%     vr = analytical.rayleigh(mat.nu,sqrt(mat.G/mat.rho)) ;
%     wr = ks.*vr ;
%     plot3(real(ks)+imag(ks),real(wr),imag(wr),':','displayname','Rayleigh') ;
% Rescaling
obj = allchild(gca) ;
[xdata,ydata] = deal(get(obj,'XData'),get(obj,'YData')) ;
if 1 % normalized k
    set(obj,{'xdata'},cellfun(@(x)x./kc,xdata,'uni',false))
    set(obj,{'ydata'},cellfun(@(y)y./2/pi/fc,ydata,'uni',false))
    xlabel 'Normalized Wavenumber $\frac{k \times h}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2h}{c_t}$'
    set(gca,'xlim',[-.25 1],'ylim',[0 1.1]) ;
else % phase velocity c
    ct = sqrt(mat.G./mat.rho) ;
    cR = analytical.rayleigh(mat.nu) ;
    cE = analytical.edge(mat.nu) ;
    obj = allchild(gca) ;
    set(obj,{'xdata'},cellfun(@(x,y)y./2/pi./fc,xdata,ydata,'uni',false))
    set(obj,{'ydata'},cellfun(@(x,y)y./x./ct,xdata,ydata,'uni',false))
%     xlabel 'Normalized Wavenumber $\frac{k \times h}{\pi}$' ; 
%     ylabel 'Normalized Frequency $\frac{f \times 2h}{c_t}$'
    set(gca,'xlim',[0 5],'ylim',[0.01 2]) ;
    plot(get(gca,'xlim'),[1 1]*cR,':') ;
    plot(get(gca,'xlim'),[1 1]*cE,':') ;
end
    
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





