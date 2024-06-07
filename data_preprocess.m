% CLEAN DATA

[file,path] = uigetfile('*.mat','SELECT A RAW DATA FILE') ;
if path==0 ; return ; end

DATA = load([path filesep file]) ;
DATA = DATA.FileData ;

%% Process

% Reshape as 2D meshgrid
nX = find(DATA.X(2:end)==DATA.X(1),1,'first') ;
nY = numel(DATA.X)/nX ;
X = reshape(DATA.X,[nX nY]) ;
Y = reshape(DATA.Y,[nX nY]) ;
f = DATA.corrFreq ;
S = reshape(DATA.AvgH1dZ,nX,nY,[]) ;


%% SAVE
[file_s,path_s] = uiputfile('*.mat','SAVE THE CLEANED DATA FILE',[path filesep file]) ;
if path_s==0 ; return ; end

save([path_s,filesep,file_s] ...
        ,'S','X','Y','f' ...
    ) ;



%% TEST
clc,clear all

[file,path] = uigetfile('*.mat','SELECT A DATA FILE') ;
if path==0 ; return ; end
load([path filesep file]) ;

%% Vectorized Wavenumber Analysis
indF = 2:10:find(f<=600e3,1,'last') ;
fun = 'cos' ;
R = 15 ;
NFFT = 1001 ;
dx = abs(X(2,1)-X(1,1)) ;

% Estimate the subspaces
St = S(:,:,indF) ;
[W,lmbda] = hrwa.subspace(St,R,fun) ;

% MUSIC pseudo-spectrum
Fw = hrwa.music(W,NFFT,fun) ;

% ESPRIT wavenumbers
ke = hrwa.esprit(W,fun) ;
ke = ke/dx ;

% Wavenumber spectrum
switch fun
    case 'exp'
        k = (1-1/NFFT)*pi/dx*linspace(-1,1,NFFT) ;
        Fw = fftshift(Fw,1) ;
    case 'cos'
        k = (1-1/NFFT)*pi/dx*linspace(0,1,NFFT) ;
end

% Display
    clf ; axis tight
    imagesc(k,f(indF)/1000,log10(abs(Fw(:,:))')) ; 
    colormap gray
    caxis([-3 0])
    plot3(real(ke),f(indF)/1e3,imag(ke),'.r','markersize',3)
    xlabel('Wavenumber $k$ (rad/mm)')
    ylabel('Frequency $\omega$ (kHz)')
    
%% REGULARIZED VERSION
indF = find(f<=200e3,1,'last'):100:find(f<=300e3,1,'last') ;
fun = 'cos' ;
R = 100 ;
NFFT = 1001 ;
beta = 100 ;

% Estimate the subspaces
St = S(:,:,indF) ; 
[q,p,nF] = size(St,1:3) ;
St = permute(St(:,:,:),[1 4 3 2]) ; % [q 1 nF p]
Fst = fft(St,[],1) ; % [q 1 nF p]
% U [m k nF]
switch fun
    case 'exp'
        m = floor(q/2) ;
        opCss = @(U)math.expcovtimes(Fst,U) ;
    case 'cos'
        m = floor(q/3) ; 
        opCss = @(U)math.coscovtimes(Fst,U) ;
end
% Regularization
    if 1 % Difference matrix
        D = eye(m*nF) - circshift(eye(m*nF),m,2) ;
        D(end-m+1:end,:) = [] ; % no circular difference
        D2 = beta*(D'*D) ;
    else % Correlation matrix
        D = diag(ones(m*(nF-1),1),m) + diag(ones(m*(nF-1),1),-m) ;
        D2 = beta*D ;
    end
% Use the Lanczos algo
    n_k_p = [m*nF R 1] ;
    opA = @(U)muticovtimesreg(U,m,opCss,0*D2) ;
    %[W,lmbda] = math.eigh(opA,n_k_p,1e-6) ;
    [W,lmbda] = eigs(opA,m*nF,eye(m*nF)+1*D2,R,'lr') ;
    lmbda = real(diag(lmbda)) ;
    W = permute(reshape(W,[m nF R]),[1 3 2]) ;
    clf ; plot(real(lmbda),'o-') ;

% MUSIC pseudo-spectrum
Fw = hrwa.music(W,NFFT,fun) ;
% ESPRIT wavenumbers
% ke = hrwa.esprit(W,fun) ;
% ke = ke/dx ;

% Wavenumber spectrum
dx = abs(X(2,1)-X(1,1)) ;
switch fun
    case 'exp'
        k = (1-1/NFFT)*pi/dx*linspace(-1,1,NFFT) ;
        Fw = fftshift(Fw,1) ;
    case 'cos'
        k = (1-1/NFFT)*pi/dx*linspace(0,1,NFFT) ;
end

% Display
    clf ; axis tight
    imagesc(k,f(indF)/1000,log10(abs(Fw(:,:))')) ; 
    colormap gray
    %caxis([-3 0])
%     plot3(real(ke),f(indF)/1e3,imag(ke),'.r','markersize',3)
    xlabel('Wavenumber $k$ (rad/mm)')
    ylabel('Frequency $\omega$ (kHz)')

    
function U = muticovtimesreg(U,m,opCss,D2)
% input/output U[m*nF k]
k = size(U,2) ;
nF = size(U,1)/m ;
U = reshape(U,[m nF k]) ; % [m nF k]
U = permute(U,[1 3 2]) ; % [m k nF]
U = opCss(U) ; % [m k nF]
U = permute(U,[1 3 2]) ; % [m nF k]
U = reshape(U,[m*nF k]) ; % [m*nF k]
if nargin>3 && ~isempty(D2)
    U = U - D2*U ; % [m*nF k]
end
end





