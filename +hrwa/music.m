function Ms = music(S,R,NFFT) 
% COMPUTES THE MUSIC SPECTRUM

nX = size(S,1) ;

% Hankel-Toeplitz block matrix
n = floor(nX/2) ; m = nX-n+1 ;
ii = (1:n)'+(0:m-1) ; % hankel
X = reshape(S(ii,:),n,[]) ;

% Signal subspace
[W,lmbda] = eig(X*X','vector') ;
[~,is] = sort(lmbda,'ascend') ;
W = W(:,is) ;
Fw = fft(W,NFFT,1) ; %dct(W,NFFT,1,'type',2) ;
Ms = 1./sum(abs(Fw(:,1:end-R)).^2,2) ;