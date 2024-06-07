%% WAVENUMBER SPECTRUM OF A PRISMATIC WAVEGUIDE
% Comparison of the SAFE method with k or w sweep
clc
clear all

% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*3*sqrt(5) ; 100/10 ; % width (mm)
    dx = min([geo.b geo.h])/4 ; % element size
% Material
    mat = [] ;
    mat.E = 70e3*(1+0.00i) ; % Young modulus (MPa)
    mat.nu = 30/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
% Eigenvalue options
    nModes = 50 ;
    dp = .001 ;
% Adimensionalized numbers
    fc = .5*sqrt(real(mat.G)/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Wavenumber domain
    kr = kc*(0:dp:1.1) ; % real wavenumbers
    ki = -1i*kc*(0:dp:.5) ; % imaginary wavenumbers
% Frequency domain
    wr = 2*pi*fc*(0:dp:1.1) ; % real wavenumbers
    wi = 2i*pi*fc*(0:dp:.5) ; % imaginary wavenumbers
% Build the mesh
    mesh = safe.mesh.gridmesh([geo.b geo.h],dx) ; % quad mesh
    mesh = safe.mesh.quad8(mesh) ; % serendipity element
    mesh = safe.mesh.quadrature(mesh,'quad2') ; % generate a quadrature rule
    clf ; axis equal ; safe.mesh.plotmesh(mesh) ; 
% Wavenumber sweep
    k = [flip(ki) kr] ;
    [Uk,Wk] = safe.computeW(mesh,mat,k,nModes,false) ;
    Kk = repmat(k,[nModes 1]) ;
% Frequency sweep
    w = [flip(wi) wr] ;
    [Uw,Kw] = safe.computeK(mesh,mat,w,nModes,false) ;
    Ww = repmat(w,[nModes 1]) ;
%% Sort branches
profile on
    edgK = sorting.branches(Uk,Wk) ;
    edgW = sorting.branches(Uw,Kw) ;
profile viewer
    clf ; axis equal 
    patch('vertices',[real(Kk(:))/kc imag(Kk(:))/kc real(Wk(:))/2/pi/fc],'faces',edgK,'EdgeColor','r')
    patch('vertices',[real(Kw(:))/kc imag(Kw(:))/kc real(Ww(:))/2/pi/fc],'faces',edgW,'EdgeColor','b')
    set(gca,'xlim',[0 1],'ylim',[-.5 0],'zlim',[0 1],'clim',[0 .5]) ;
    view([40 30])
    
%% Display
    clf reset ; hold on
    axis equal
    scatter3(real(Kk(:))/kc,imag(Kk(:))/kc,real(Wk(:))/2/pi/fc,10,imag(Wk(:))/2/pi/fc,'displayname','$\omega(k)$') ;
    scatter3(real(Kw(:))/kc,imag(Kw(:))/kc,real(Ww(:))/2/pi/fc,5,imag(Ww(:))/2/pi/fc,'filled','displayname','$k(\omega)$') ;
% Legends
%     xlabel 'Wavenumber $k$ (rad/mm)' ; 
%     ylabel 'Frequency $\omega$ (rad/s)'
    lgd = legend ; lgd.Location = 'southeast' ;
    set(gca,'xlim',[0 1],'ylim',[-.5 0],'zlim',[0 1],'clim',[0 .5]) ;
    view([40 30])
    colormap(interp1([0 0 1;[1 1 1]*.75;1 0 0],linspace(1,3,100)'))
    colormap(interp1([[1 1 1]*.5;1 0 0],linspace(1,2,100)'))
    colorbar
    
%% Rulers
    hAxis = gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.XRuler.SecondCrossoverValue  = 0; % X crossover with Z axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.YRuler.SecondCrossoverValue  = 0; % Y crossover with Z axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis



