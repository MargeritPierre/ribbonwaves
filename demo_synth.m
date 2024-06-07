%% SYNTHESIZE A MEASUREMENT ON A RIBBON
clc
clear all

% Ribbon geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = geo.h*sqrt(3*5) ; % width (mm)
    geo.L = 100 ; % length (mm)
    dx = geo.b/10 ; % element size
% Material
    mat = [] ;
    mat.E = 70e3*(1+0.0000i) ; % Young modulus (MPa)
    mat.nu = 30/100 ; % Poisson ratio
    mat.rho = 2700e-12 ; % material density (tons/mm^3)
    mat.xi = pi^2/12 ; % shear correction coefficient
    mat = material.coefficients(mat) ;
% Adimensionalized numbers
    fc = .5*sqrt(real(mat.G)/mat.rho)/geo.h ; % height-cutoff frequency
    kc = pi/geo.h ; % height-cutoff wavenumber
% Frequencies
    freq = linspace(0.001,1,100)*fc ;
% Build the mesh
    mesh = safe.mesh.gridmesh([geo.L geo.b],dx) ; % quad mesh
    mesh = safe.mesh.quad8(mesh) ; % serendipity element
    mesh = safe.mesh.quadrature(mesh,'quad2') ; % generate a quadrature rule
    clf ; axis equal ; safe.mesh.plotmesh(mesh) ; 
% Problem matrices
    % Degrees of freedom: u = [u1;u2;u3;phi1;phi2]
    % Build the interpolation matrices
        mesh = safe.mesh.interpolation(mesh) ; % generate interpolation matrices
        [N,Dx,Wq,O] = deal(mesh.N,mesh.Dx,mesh.Wq,mesh.O) ;
        NWN = N'*Wq*N ;
        Q = mat.Q.*[1 mat.nu 0 ; mat.nu 1 0 ; 0 0 .5*(1-mat.nu)] ;
    % In-Plane motion
        % strains Em = [E11;E22;2E12] = Bm*u
            Bm = [  ...
                    Dx{1} O O O O ; .... E11 = u1,1
                    O Dx{2} O O O ; ... E22 = u2,2
                    Dx{2} Dx{1} O O O ; ... 2E12 = u1,2+u2,1
                 ] ;
        % stiffness
            Km = geo.h*Bm'*kron(sparse(Q),Wq)*Bm ;
        % inertia
            Mm = geo.h*mat.rho*kron(sparse(diag([1 1 0 0 0])),NWN) ;
    % Bending motion
        % strains Eb = [K11;K22;2K12] = Bb*u
            Bb = [  ...
                    O O O Dx{1} O ; .... K11 = phi1,1
                    O O O O Dx{2} ; ... K22 = phi2,2
                    O O O Dx{2} Dx{1} ; ... KE12 = phi1,2+phi2,1
                 ] ;
        % stiffness
            Kb = geo.h^3/12*Bb'*kron(sparse(Q),Wq)*Bb ;
        % inertia
            Mb = geo.h^3/12*mat.rho*kron(sparse(diag([0 0 0 1 1])),NWN) ;
    % Out-of-plane shear motion
        % strains Es = [d1;d2] = Bs*u
            Bs = [  ...
                    O O Dx{1} N O ; .... d1 = phi1+u3,1
                    O O Dx{2} O N ; ... d2 = phi2+u3,2
                 ] ;
        % stiffness F
            Ks = geo.h*mat.G*mat.xi*(Bs'*kron(speye(2),Wq)*Bs) ;
        % inertia
            Ms = geo.h*mat.rho*kron(sparse(diag([0 0 1 0 0])),NWN) ;
    % Total system
        K = Km+Kb+Ks ;
        M = Mm+Mb+Ms ; M = .5*(M+M') ;
% Boundary conditions
    keepDOF = true(size(mesh.X,1),5) ;
    keepDOF(mesh.X(:,1)<=eps,:) = false ; % clamp the left side
    %keepDOF(mesh.X(:,1)>geo.L-eps,:) = false ;
    K = K(keepDOF,keepDOF) ;
    M = M(keepDOF,keepDOF) ;
% Applied load
    f = zeros(size(mesh.X,1),5) ;
    f(mesh.X(:,1)>=geo.L-eps & mesh.X(:,2)<=eps,1:3) = 1 ; % 3D load in the right-bottom corner
    f = f(keepDOF) ;
% Solve for all frequencies
    U = zeros(size(mesh.X,1)*5,numel(freq)) ;
    for fr = 1:numel(freq)
        U(keepDOF,fr) = (K-(2*pi*freq(fr))^2*M)\f ;
        fr
    end
    U = reshape(U,size(mesh.X,1),5,numel(freq)) ;
%% Display
    clf ; axis equal
    p = safe.mesh.plotmesh(mesh,100*real(U(:,1:3,end))) ;
    p.Marker = 'none' ; p.EdgeColor = 'none' ;
%% APPLY HRWA
    indX = 1:prod(mesh.nX+1) ;
    R0 = 1:20 ; m = 6 ;
    k = NaN(numel(freq),max(R0)) ;
    dx = geo.L/mesh.nX(1) ;
    for fr = 1:numel(freq)
        u = U(indX,1:3,fr) ;
        u = reshape(u,mesh.nX(2)+1,mesh.nX(1)+1,[]) ;
        u = u(:,1+m:end-m,:) ;
        out = ESPRIT_fcn(u...
                            ,'R0',R0 ...
                            ,'DIMS_K',2 ...
                            ,'FUNC','cos' ...
                            ,'CRITERION','ESTER' ...
                            ,'M/L',.5 ...
                            ,'DECIM',[1 1 1] ...
                            ) ;
        k(fr,1:numel(out.K)) = out.K/dx ;
    end
    clf ; plot3(real(k),freq(:),imag(k),'.') ;
    
%% Ribbon Model determinant
    nPix = 500 ;
    [D,Kd,Wd,colors] = ribbon.detMap(geo,mat,[-.25i 0 1]*kc,2*pi*fc*[0 1.1],nPix) ;
% Display
    clf ; hold on
    surf(real(Kd)+0*imag(Kd),Wd,D(:,:,1)*0,D,'facecolor','interp','edgecolor','none','handlevisibility','off') ;
    plot(NaN,NaN,'color',colors(1,:),'displayname','IP-S') ;
    plot(NaN,NaN,'color',colors(2,:),'displayname','IP-A') ;
    plot(NaN,NaN,'color',colors(3,:),'displayname','OP-S') ;
    plot(NaN,NaN,'color',colors(4,:),'displayname','OP-A') ;
    plot(real(k)+0*imag(k),2*pi*freq(:),'.') ;
        
         