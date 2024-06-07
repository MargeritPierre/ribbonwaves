%% PARAMETRIC ESTIMATION OF A COMBINATION OF COSINE FUNCTIONS
% Signal model: s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
clc
clear all

% Parameters
R = 10 ; % signal order
nX = 300 ; % signal length along wavenumber dimension
nY = 6 ; % number of points on the "mode"
nF = 5 ; % Number of "frequencies"
SNR = 1e10000 ; % signal-to-noise ratio

% Wavenumbers, amplitudes & phases
k = .75*pi*rand(R,nF).*(1-0.001i*rand(R,nF)) ; 2*pi.*linspace(.5/R,1-.5/R,R)' ;
A = randn(R,nY,nF) + 1i*randn(R,nY,nF) ;
PHI = pi/2*(2*rand(R,nY,nF)-1) ;
% clf ; axis equal ; plot(exp(1i*k),'o')

% Signal generation; size(S) = [nX nY]
x = (0:nX-1)' ;
S = sum(permute(A,[4 2 3 1]).*exp(-1i*permute(k,[3 4 2 1]).*x+permute(PHI,[4 2 3 1])),4) ; % pure signal
B = randn(nX,nY) + 1i*randn(nX,nY) ; % noise
S = S + B*norm(S(:))/norm(B(:))/SNR ;
% clf ; plot(real(S(:,:))) ; plot(imag(S(:,:))) ;

%% Estimate the signal subspace
% profile on 
tic ; [W,lmbda] = hrwa.subspace(S,R) ; toc
% profile viewer
tic
Wref = NaN(size(W)) ;
lmbdaref = NaN(size(lmbda)) ;
for ff = 1:nF
    H = hankel(1:size(W,1),size(W,1):size(S,1)) ;
    H = reshape(S(H,:,ff),size(H).*size(S(1,:,ff))) ;
    [Wref(:,:,ff),lm] = eigs(H*H',R,'lm') ;
    lmbdaref(:,:,ff) = diag(lm) ;
end
toc

%% Compute the MUSIC spectrum
NFFT = 10*nX ;
Ms = hrwa.music(S,2*R,NFFT) ;

% Apply ESPRIT
[ke,Ae,PHIe] = hrwa.esprit(S,R) ;

% Estimation errors
% Sort estimations
errMatK = abs(ke(:)-k(:).') ;
[~,cloK] = min(errMatK,[],1) ; % find the closest estimation
ke = ke(cloK) ;
Ae = Ae(cloK,:) ;
PHIe = PHIe(cloK,:) ;
% errors
errK = norm(ke-k) 
errA = norm(Ae-A)
errPHI = norm(PHIe-PHI)

% Display
clf ; 
plot((0:NFFT-1)/NFFT*2*pi,log10(Ms))
plot(real(k(:))'.*[1;1],get(gca,'ylim')'.*ones(1,R),'-','linewidth',.1)

