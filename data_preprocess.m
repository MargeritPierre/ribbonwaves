%% CLEAN RAW DATA

[file,path] = uigetfile('*.mat','SELECT A RAW DATA FILE') ;
if path==0 ; return ; end

DATA = load([path filesep file]) ;
DATA = DATA.FileData ;

% Reshape as 2D meshgrid
f = DATA.corrFreq ;
isXscan = DATA.X(2)~=DATA.X(1) ;
if isXscan
    nX = find(DATA.X(2:end)==DATA.X(1),1,'first') ;
    nY = numel(DATA.X)/nX ;
    sz = [nX nY] ;
else
    nY = find(DATA.Y(2:end)==DATA.Y(1),1,'first') ;
    nX = numel(DATA.X)/nY ;
    sz = [nY nX] ; 
end
X = reshape(DATA.X,sz) ;
Y = reshape(DATA.Y,sz) ;
S = reshape(DATA.AvgH1dZ,[sz numel(f)]) ;

% Rotate so that the longest direction is in first dimension
isXlongest = range(X(:))>range(Y(:)) ;
if isXlongest && isXscan % ok
elseif isXlongest && ~isXscan % permute X and Y
    X = permute(X,[2 1 3]) ;
    Y = permute(Y,[2 1 3]) ;
    S = permute(S,[2 1 3]) ;
elseif ~isXlongest && isXscan % swap X and Y then permute
    [X,Y] = deal(Y,X) ;
    X = permute(X,[2 1 3]) ;
    Y = permute(Y,[2 1 3]) ;
    S = permute(S,[2 1 3]) ;
elseif ~isXlongest && ~isXscan %  swap X and Y
    [X,Y] = deal(Y,X) ;
end

% SAVE
[file_s,path_s] = uiputfile('*.mat','SAVE THE CLEANED DATA FILE',[path filesep file]) ;
if path_s==0 ; return ; end

save([path_s,filesep,file_s] ...
        ,'S','X','Y','f' ...
        ,'-v7.3'...
    ) ;



%% LOAD CLEAN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear all

[file,path] = uigetfile('*.mat','SELECT A DATA FILE') ;
if path==0 ; return ; end
load([path filesep file]) ;

%% Vectorized Wavenumber Analysis
fmin = 100 ; fmax = 200e3 ; nF = 1000 ;
indF = unique(round(linspace(find(f>=fmin,1,'first'),find(f<=fmax,1,'last'),nF))) ;
fun = 'cos' ;
R = 20 ;
NFFT = 1001 ;

% Select signal samples
disp('--- SIGNAL TRUNCATION ---') ;
St = S(:,:,indF) ;

% Estimate the subspaces
disp('--- SUBSPACE ESTIMATION ---') ;
[W,lmbda] = hrwa.subspace(St,R,fun) ;
W = W.*permute(sqrt(lmbda),[2 1 3]) ; % keep the eigenvalue info in the signal subspace

% MUSIC pseudo-spectrum
disp('--- MUSIC SPECTRUM ---') ;
Fw = hrwa.music(W,NFFT,fun) ;

% Wavenumber spectrum
dx = abs(X(2,1)-X(1,1)) ;
switch fun
    case 'exp'
        k = (1-1/NFFT)*pi/dx*linspace(-1,1,NFFT) ;
        Fw = fftshift(Fw,1) ;
    case 'cos'
        k = (1-1/NFFT)*pi/dx*linspace(0,1,NFFT) ;
end

% ESPRIT wavenumbers
    disp('--- ESPRIT WAVENUMBERS ---') ;
    [ke,Ue,Sre] = hrwa.esprit(W,fun,St) ;
    ke = ke/dx ;
    Se = sum(Sre,4) ; % approximation of St
    Ee = permute(sum(abs(Sre).^2,1:2),[4 3 1 2]) ; % component energies [R nF]
    if strcmp(fun,'cos') ; Ee = Ee(1:R,:) + Ee(R+(1:R),:) ; end
    Eer = Ee./sum(Ee,1) ; % relative energies [R nF]
    fe = repmat(f(indF),[R 1]) ; % frequencies corresponding to the estimated wavenumbers


% DISPLAY
    clf reset ; axis tight
    xlabel('Wavenumber $Re(k)$ (rad/mm)')
    ylabel('Frequency $\omega$ (kHz)')
    zlabel('Wavenumber $Im(k)$ (rad/mm)')
    set(gca,'zdir','reverse') ;
    %set(gca,'xlim',[0 2]);
    % MUSIC
        clim = [-4 0] ;
        ms = log10(abs(Fw(:,:)))' ;
        ms = interp1(clim,[0 1],max(clim(1),min(ms,clim(2)))) ;
        ms = repmat(ms,[1 1 3]) ;
        image(k,f(indF)/1e3,ms) ; 
    % ESPRIT
        markersize = 50 ;
        cdata = Eer ; % color with the relative energy of the signal components
        sc = scatter3(real(ke(:)),fe(:)/1e3,imag(ke(:)),markersize,cdata(:),'.') ;
        set(gca,'colorscale','log') ;
        ylabel(colorbar('location','east'),'Relative Energy','interpreter','latex') ;

%% Clean ESPRIT RESULTS    
    minRelativeEnergy = 1e-2 ; % minimum relative energy in the signal
    maxImagK = 1e0*pi/dx ; % maximum imaginary part of the wavenumbers
    maxSpatialDecay = 1e0 ; % maximum allowed spatial decay
    remPositiveDecay = true ; % remove wavenumbers with positive decay
    
    sc.XData = real(ke(:)) ;
    valid = true(size(ke)) ;
    valid = valid & Eer>minRelativeEnergy ;
    valid = valid & abs(imag(ke))<maxImagK ;
    valid = valid & abs(imag(ke)./real(ke))<maxSpatialDecay ;
    if remPositiveDecay ; valid = valid & imag(ke)./real(ke)<=0 ; end
    sc.XData(~valid) = NaN ;
    
    Kv = ke(valid) ; Fv = fe(valid) ; Ev = Eer(valid) ;
    
%% PHASE VELOCITY
    clf ; 
    plot3(Fv,2*pi*Fv./real(Kv)/1000,2*pi*Fv./imag(Kv)/1000,'.k') ;
    %set(gca,'xscale','log','yscale','log')
    xlabel 'Frequency $\omega$ (Hz)' ;
    ylabel 'Phase Velocity $c$ (m/s)'
    
    
%% INITIAL COMPARISON WITH THE RIBBON MODEL
% Geometry
    geo = [] ;
    geo.h = 20/10 ; % height (mm)
    geo.b = 195/10 ; % width (mm)
% Material
    mat = [] ;
    mat.E = 52e2*(1+0.04i) ; % Young modulus (MPa)
    mat.nu = 42/100 ; % Poisson ratio
    mat.rho = 1200e-12 ; % material density (tons/mm^3)
    mat = material.coefficients(mat) ;
% Initial solutions
    wm = 2*pi*f(indF(1:10:end)) ;
    [Kips,Kipa,Kops,Kopa] = ribbon.computeK(geo,mat,wm) ;
    Wm = repmat(wm,[size(Kips,1) 1]) ;
% Display
    clf ;
    plot3(real(Kips(:)),Wm(:)/2/pi/1000,imag(Kips(:)),'displayname','$k_{IP}^{S}$','tag','model') ;
    plot3(real(Kipa(:)),Wm(:)/2/pi/1000,imag(Kipa(:)),'displayname','$k_{IP}^{A}$','tag','model') ;
    plot3(real(Kops(:)),Wm(:)/2/pi/1000,imag(Kops(:)),'displayname','$k_{OP}^{S}$','tag','model') ;
    plot3(real(Kopa(:)),Wm(:)/2/pi/1000,imag(Kopa(:)),'displayname','$k_{OP}^{A}$','tag','model') ;
    set(findobj(gca,'tag','model'),'linestyle','none','marker','o','linewidth',1,'markersize',5) ;
    plot3(real(Kv),Fv/1000,imag(Kv),'.k','displayname','HRWA','markersize',5) ;
    set(gca,'zlim',[-1 1]*.2*pi/geo.b,'zdir','reverse')
    legend('location','southeast')
    
    
    
%% STOP HERE: WHAT FOLLOWS IS EXPERIMENTAL :)
    
    
    
    
    
    
%% Find the closest roots from the HRWA results
    Dv = ribbondet(geo,mat,Kv,2*pi*Fv) ;
    Kd = math.muller(@(K)ribbondet(geo,mat,K,2*pi*Fv,false),Kv,100) ;
    %Kd = math.gradientdescent(@(K)ribbondet(geo,mat,K,2*pi*Fv),Kv,100) ;
    Dd = ribbondet(geo,mat,Kd,2*pi*Fv) ;
    vk = [Kv Kd Kd*NaN].' ; vf = [Fv Fv Fv*NaN]' ;
    
    clf ;
    plot3(real(vk(:)),vf(:)/1000,imag(vk(:)),'b','displayname','link','linewidth',.1) ;
    plot3(real(Kv),Fv/1000,imag(Kv),'.k','displayname','HRWA') ;
    plot3(real(Kd),Fv/1000,imag(Kd),'.r','displayname','closest') ;
    legend('location','southeast')
    %%
    clf
    nPix = 500 ;
    kk = linspace(0,max(abs(Kv(:))),nPix) ;
    ff = linspace(0,max(abs(Fv(:))),nPix) ; 
    tic ; dd = ribbondet(geo,mat,kk+0*ff(:),2*pi*ff(:)+kk,true) ; toc
    dp = dd ; 
    %dp = log10(abs(dd)+1).*exp(1i*angle(dd)) ;
    dp = real(dp) ;
    surf(kk+0*ff(:),ff(:)+kk,dp,'edgecolor','none','facecolor','interp') ;
    %set(gca,'colorscale','log')
    
    
%% Util functions
    function D = ribbondet(geo,mat,k,w,normalize)
    % Full determinant of the ribbon model
        delta = 1e-6 ; n = 3 ; 
        if nargin<5 ; normalize = false ; end
        if normalize % normalize the determinant with the local values
        % size of the local zone
            dk = delta*max(abs(k(:)))*linspace(-1,1,n) ;
            dw = delta*max(abs(w(:)))*linspace(-1,1,n) ;
            [~,DK,DW] = ndgrid(1,dk,dw) ;
        % Match sizes (so that size(k)==size(w))
            k = k + 0*w ; 
            w = w + 0*k ;
        % Add the local zone
            K = k(:) + DK ;
            W = w(:) + DW ;
        else
            K = k ; W = w ;
        end
    % Compute the determinant
        if 0
            m = ribbon.model(geo,mat,K,W) ;
            D = m.Dips_n.*m.Dipa_n.*m.Dops_n.*m.Dopa_n ;
        else
            vt = sqrt(mat.G./mat.rho) ;
            [~,~,~,~,D] = ribbon.normdet(geo,mat,.5.*geo.b.*K,.5.*geo.b.*W./vt) ;
        end
    % Normalize if needed
        if normalize
            D = D(:,:) ;
            %D = mean(D,2)./mean(abs(D),2) ;
            D = D(:,ceil(end/2))./mean(abs(D),2) ;
            D = reshape(D,size(k)) ;
        end
    end




