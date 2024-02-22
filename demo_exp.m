%% PARAMETRIC ESTIMATION OF A COMBINATION OF COSINE FUNCTIONS
% Signal model: s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
clc
clear all

% Parameters
R = 1 ; % signal order
nX = 100 ; % signal length along wavenumber dimension
nY = 1 ; % number of points on the "mode"
SNR = 1e1 ; % signal-to-noise ratio

% Wavenumbers, amplitudes & phases
k = .75*pi*rand(R,1).*(1-0.05i*rand(R,1)) ;
A = randn(R,nY) + 1i*randn(R,nY) ;
PHI = pi/2*(2*rand(R,nY)-1) ; pi*rand(R,nY) ;
clf ; axis equal ; plot(exp(1i*k),'o')

% Signal generation; size(S) = [nX nY]
x = (0:nX-1)' ;
S = sum(permute(A,[3 2 1]).*cos(permute(k,[3 2 1]).*x+permute(PHI,[3 2 1])),3) ; % pure signal
B = randn(nX,nY) + 1i*randn(nX,nY) ; % noise
S = S + B*norm(S(:))/norm(B(:))/SNR ;
clf ; plot(real(S)) ; plot(imag(S)) ;

% Compute the MUSIC spectrum
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

