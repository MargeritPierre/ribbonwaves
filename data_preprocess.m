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

%% DISPLAY DEFORMED SHAPES

ff = 203 ; 
phi = 352*pi/18 ;
amp = .5*norm([range(X(:)) range(Y(:))]) ;

uu = S(:,:,ff) ;
uu = amp*uu./max(abs(uu(:))) ;
uu = real(uu.*exp(1i*phi)) ;

clf ; axis equal tight
surf(X,Y,uu,'facecolor','interp') ;
view([-20 80])
set(gca,'zlim',amp*[-1 1])
title(num2str(f(ff)/1000,5)+" kHz")


%% Vectorized Wavenumber Analysis
fmin = 0.01e3 ; fmax = 70e3 ; nF = 1000 ;
indF = unique(round(linspace(find(f>=fmin,1,'first'),find(f<=fmax,1,'last'),nF))) ;
fun = 'cos' ;
R = 3 ;
NFFT = 1001 ;
nSmooth = 0 ; smoothing = 'cat' ;

% Select signal samples
disp('--- SIGNAL TRUNCATION ---') ;
ss = -nSmooth:nSmooth ;
St = reshape(S(:,:,indF+ss(:)),size(S,1),size(S,2),numel(ss),numel(indF)) ;
switch smoothing
    case 'cat' % concatenate neightboring frequencies as more snapshots
        St = reshape(St,size(S,1),[],numel(indF)) ;
    case 'corr' % compute the correlation
        nCorr =(numel(ss)*size(S,1))-1 ; 
        if nSmooth>0
            %St = St.*blackman(size(S,1)) ;
            St = fft(St,nCorr,1) ;
            St = prod(St,3) ;
            St = ifft(St,[],1) ;
        end
        St = St(:,:,:) ;
end

% Estimate the subspaces
disp('--- SUBSPACE ESTIMATION ---') ;
[W,lmbda] = hrwa.subspace(St,R,fun) ;
%W = W.*permute(sqrt(lmbda),[2 1 3]) ; % keep the eigenvalue info in the signal subspace

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
    %cla
    xlabel('Wavenumber $Re(k)$ (rad/mm)')
    ylabel('Frequency $\omega$ (kHz)')
    zlabel('Wavenumber $Im(k)$ (rad/mm)')
    set(gca,'zdir','reverse') ;
    %set(gca,'zlim',[-.02 0])
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
    minRelativeEnergy = 20e-3 ; % minimum relative energy in the signal
    maxImagK = 1e0*pi/dx ; % maximum imaginary part of the wavenumbers
    maxSpatialDecay = 20e-2 ; % maximum allowed spatial decay
    remPositiveDecay = true ; % remove wavenumbers with positive decay
    
    sc.XData = real(ke(:)) ;
    valid = true(size(ke)) ;
    valid = valid & Eer>minRelativeEnergy ;
    valid = valid & abs(imag(ke))<maxImagK ;
    valid = valid & abs(imag(ke)./real(ke))<maxSpatialDecay ;
    if remPositiveDecay ; valid = valid & imag(ke)./real(ke)<=0 ; end
    sc.XData(~valid) = NaN ;
    
    Kv = ke(valid) ; Fv = fe(valid) ; Ev = Eer(valid) ;
    
%% INCREASE ACCURACY WITH LSFIT
    Kf = ke ;
    for ff = 1:1:numel(indF) 
        [Kf(valid(:,ff),ff),Uf] = hrwa.lsfit(St(:,:,ff),Kf(valid(:,ff),ff)*dx,fun) ;
    end
    Kf = Kf(valid)/dx ;
    delete(findall(gcf,'tag','LSFIT'))
    plot3(real(Kf),Fv/1e3,imag(Kf),'.b','tag','LSFIT') ;
    set(gca,'zlim',[-0.1 0])
%% PHASE VELOCITY
    clf ; 
    plot3(Fv,2*pi*Fv./real(Kv)/1000,2*pi*Fv./imag(Kv)/1000,'.k') ;
    %set(gca,'xscale','log','yscale','log')
    xlabel 'Frequency $\omega$ (Hz)' ;
    ylabel 'Phase Velocity $c$ (m/s)'
    
    
%% INITIAL COMPARISON WITH THE RIBBON MODEL
% Geometry
    geo = [] ;
    geo.h = 75/100 ; % height (mm)
    geo.b = 40/10 ; % width (mm)
% Material
    mat = [] ;
    mat.E = 15e2*(1+0.04026i) ; % Young modulus (MPa)
    mat.nu = 40/100 ; % Poisson ratio
    %mat.G = 1212*(1+0.0438i) ; % Poisson ratio
    mat.rho = 1120e-12 ; % material density (tons/mm^3)
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
    
    
    
    
    
    
%% ESTIMATION OF THE RIBBON PARAMETERS
    % Initial arguments
        args0 = struct(...
                        ... 'h' , 2 , 'b' , 10 , 'E' , 71e3 , 'G', 27e3 , 'rho', 2700e-12 ... Alu 10x2
                        ... 'h' , 2 , 'b' , 15 , 'E' , 71e3 , 'G', 27e3 , 'rho', 2700e-12 ... Alu 15x2
                        ...'h' , 2 , 'b' , 19 , 'E' , 51e2 , 'G', 19e2 , 'rho', 1200e-12 ... % PVC 19x2
                        'h' , .75 , 'b' , 4 , 'E' , 3.8e3 , 'G', 1.43e3 , 'rho', 1100e-12 ... % VW100 4x0.75
                        ...'h' , 1 , 'b' , 8 , 'E' , 3.8e3 , 'G', 1.43e3 , 'rho', 1100e-12 ... % VW100 8x1
                        ...'h' , .75 , 'b' , 4 , 'E' , 2.6e3 , 'G', .9e3 , 'rho', 1100e-12 ... % VW66 4x0.75
                        ...'h' , .75 , 'b' , 4 , 'E' , 1.2e3 , 'G', .35e3 , 'rho', 1100e-12 ... % EB100 4x0.75
                      ) ;
        optim = {...
                    'E' ...
                    , 'G' ...
                    ..., 'h' ...
                    ..., 'b' ...
                } ;
        conFun = @(args)args ; @(args)struct(... impose constraints
                            'h',real(args.h) ...
                            ,'b',real(args.b) ...
                            ,'E',args.E+1i*abs(imag(args.E)) ...
                            ,'G',real(args.G)+1i*abs(imag(args.G)) ...
                            ,'rho',real(args.rho) ...
                            ) ;
        args = conFun(args0) ;
    % Descent algorithm parameters
        normalize = true ; % normalized determinant ?
        delta = 1e-6 ; % for partial derivative computation
        sigmaK = 5e-2*pi./args0.h ; % confidence interval (data weighting)
    % Functions
        geoFun = @(args)struct('h',args.h,'b',args.b) ;
        matFun = @(args)material.coefficients(struct('E',args.E,'G',args.G,'rho',args.rho)) ;
        detFun = @(args,K,W)ribbondet(geoFun(args),matFun(args),K,W,normalize) ;

%% Find the closest roots from the HRWA results
    Kd = math.muller(@(K)detFun(args,K,2*pi*Fv),Kv) ;
    vk = [Kv Kd Kd*NaN].' ; vf = [Fv Fv Fv*NaN]' ;
    
% Compute the Jacobian
    Dd = detFun(args,Kd,2*pi*Fv) ;
% Derivatives of the determinant w.r.t parameters
    dD_dp = zeros(numel(Dd),numel(optim)) ;
    for pp = 1:numel(optim)
        dp = abs(args0.(optim{pp})*delta) ;
        dargs = setfield(args,optim{pp},args.(optim{pp})+dp) ;
        dD_dp(:,pp) = (detFun(dargs,Kd,2*pi*Fv)-Dd)./dp ;
    end
% Derivatives of the determinant w.r.t wavenumber
    dk = delta.*max(abs(Kd(:))) ;
    dD_dk = (detFun(args,Kd+dk,2*pi*Fv)-Dd)./dk ;
% Jacobian
    Jp = -dD_dp./dD_dk ; % J = dk/dp
    
% Display the fit
    fs = .5*sqrt(real(args0.G)/args0.rho)/args0.h ; % height-cutoff frequency
    ks = pi/args0.h ; % height-cutoff wavenumber
    cla ;
    plot3(real(vk(:))./ks,vf(:)./fs,imag(vk(:))./ks,'b','displayname','link','linewidth',.1) ;
    plot3(real(Kv)./ks,Fv./fs,imag(Kv)./ks,'.k','displayname','$k_e$') ;
    plot3(real(Kd)./ks,Fv./fs,imag(Kd)./ks,'.r','displayname','$k_e^D$') ;
%     for pp = 1:numel(optim)
%         quiver3( ...
%                     real(Kd)./ks,Fv./fs,imag(Kd)./ks ...
%                     ,real(Jp(:,pp))./ks,0*Fv./fs,imag(Jp(:,pp))./ks ...
%                     ,'displayname',"$dD/d"+optim{pp}+"$" ...
%                     ,'AutoScaleFactor',10 ...
%                 ) ;
%     end
    legend('location','southeast')
    
% Update
    w = Ev.*exp(-(abs(Kd-Kv)./sigmaK).^2) ;
    j = Jp'*(w.*(Kd-Kv)) ;
    H = Jp'*(w.*Jp) ;
    dp = -H\j ;
    for pp = 1:numel(optim)
        args.(optim{pp}) = args.(optim{pp}) + dp(pp) ;
    end
    args = conFun(args) ;
    args
   
    
%% Subspace similarity

dW = W(:,:,1:end-1)-W(:,:,2:end) ;
dW = sum(abs(dW).^2,2) ; 
dW = dW(:,:) ;

clf ; 
mysubplot(1,2,1) ;
    clim = [-2 0] ;
    ms = log10(abs(Fw(:,:)))' ;
    ms = interp1(clim,[0 1],max(clim(1),min(ms,clim(2)))) ;
    ms = repmat(ms,[1 1 3]) ;
    image(k,f(indF)/1e3,ms) ; 
mysubplot(1,2,2) ;
%     imagesc(k,f(indF)/1e3,squeeze(abs(W(:,1,:))).')
    surf(squeeze(abs(W(:,1,:))).','facecolor','interp')
    
    %%
    Ss = St ;
    Ss = Ss./sqrt(sum(abs(Ss).^2,1:2)) ;
    Ss = real(Ss) ;
    Ss = sum(Ss(:,2,:),2) ;
    clf ; surf(squeeze(Ss).','facecolor','interp')
    
%% REGULARIZED SUBSPACE 1

R = 2 ;
beta = 1e3 ;
NFFT = 1001 ;
[nS,p,nF] = size(St) ;

Fs = fft(permute(St(:,:,:),[1 4 3 2]),[],1) ; % [q 1 nF p]

switch fun
    case 'exp'
        m = floor(nS/2) ;
        opCss = @(U)math.expcovtimes(Fs,U) ;
    case 'cos'
        m = floor(nS/3) ;
        opCss = @(U)math.coscovtimes(Fs,U) ;
end
    
% Difference matrix
D = eye(m*nF) - circshift(eye(m*nF),m,2) ;
D(end-m+1:end,:) = [] ; % no circular difference

opA = @(U) reshape(permute(opCss(permute(reshape(U,m,nF,[]),[1 3 2])),[1 3 2]),m*nF,[]) ;
opB = eye(m*nF)+beta*(D'*D) ;

% Eigenspace
[Wr,s] = eigs(opA,m*nF,opB,R,'lm') ;
[Wr,RR] = qr(Wr,0) ;
[s,is] = sort(diag(s),'descend') ;
Wr = reshape(Wr(:,is),[m nF R]) ;
Wr = Wr./sqrt(sum(abs(Wr).^2,1)) ;
%clf ; surf(real(W(:,:,1)),'facecolor','interp') ;
Wr = permute(Wr,[1 3 2]) ;


% MUSIC pseudo-spectrum
Fwr = hrwa.music(Wr,NFFT,fun) ;
if strcmp(fun,'exp') ; Fwr = fftshift(Fwr,1) ; end

clf ; 
mysubplot(1,2,1) ;
    clim = [-2 0] ;
    ms = log10(abs(Fw(:,:)))' ;
    ms = interp1(clim,[0 1],max(clim(1),min(ms,clim(2)))) ;
    ms = repmat(ms,[1 1 3]) ;
    image(k,f(indF)/1e3,ms) ; 
mysubplot(1,2,2) ;
    clim = [-2 0] ;
    ms = log10(abs(Fwr(:,:)))' ;
    ms = interp1(clim,[0 1],max(clim(1),min(ms,clim(2)))) ;
    ms = repmat(ms,[1 1 3]) ;
    image(k,f(indF)/1e3,ms) ;
    
    
%% Util functions
    function D = ribbondet(geo,mat,k,w,normalize)
    % Full determinant of the ribbon model
        delta = 1e-4 ; n = 2 ; 
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
            [Dips,Dipa,Dops,Dopa,D] = ribbon.normdet(geo,mat,.5.*geo.b.*K,.5.*geo.b.*W./vt) ;
            %D = Dops.*Dopa ;
        end
    % Normalize if needed
        if normalize
            D = D(:,:) ;
            %D = mean(D,2)./mean(abs(D),2) ;
            D = D(:,ceil(end/2))./mean(abs(D),2) ;
            D = reshape(D,size(k)) ;
        end
    end




