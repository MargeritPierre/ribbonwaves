%% WAVENUMBER SPECTRUM OF A PRISMATIC WAVEGUIDE
% Comparison between the SAFE method and the reduced ribbon model
clc
clear all

% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*5 ; 100/10 ; % width (mm)
    dx = min([geo.b geo.h])/4 ; % element size
% Material
    mat = [] ;
    mat.E = 3.5e3*(1+0.0000i) ; % Young modulus (MPa)
    mat.nu = 30/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
% Adimensionalized numbers
    fc = .5*sqrt(mat.G/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Wavenumber domain
    kr = kc*linspace(0.01,1,100) ; % real wavenumbers
    ki = -1i*kc*linspace(0.01,.25,25) ; % imaginary wavenumbers
% Ribbon Model determinant
    nPix = 500 ;
    f = [0 1.1]*fc ; % frequency range
    kd = [ki(end) 0 kr(end)] ; % wavenumber line
    [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,kd,2*pi*f,nPix) ;
% Eigenvalue options
    nModes = 25 ;
% Build the mesh
    mesh = safe.mesh.gridmesh([geo.b geo.h],dx) ; % quad mesh
    clf ; axis equal ; safe.mesh.plotmesh(mesh) ; 
% Compute the SAFE frequencies & modes
    [Ur,wr] = safe.computeW(mesh,mat,kr,nModes) ;
    [Ui,wi] = safe.computeW(mesh,mat,ki,nModes) ;
    Ks = repmat([real(kr(:));imag(ki(:))],[1 nModes]) ;
    Ws = [real(wr) ; real(wi)] ;
% Display
    clf ; hold on
    surf(real(Kd)+imag(Kd),Wd,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP-S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP-A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP-S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP-A') ;
    plot(Ks(:),Ws(:),'.k','displayname','SAFE') ; 
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
    set(gca,'xlim',[-.25 1],'ylim',[0 1.1]) ;
    
    



