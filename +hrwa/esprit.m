function [k,A,PHI] = esprit(S,R)
% PARAMETRIC ESTIMATION OF A COMBINATION OF COSINE FUNCTIONS
% Signal model: s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }

nX = size(S,1) ;

% Hankel-Toeplitz block matrix
n = floor(nX/3) ; m = nX-2*n+1 ;
i1 = (n+1:2*n)'+(0:m-1) ; % hankel
i2 = flip(1:n)'+(0:m-1) ; % toeplitz
X = reshape(S(i1,:) + S(i2,:),n,[]) ;

% Signal subspace
[W,lmbda] = eig(X*X','vector') ;
[lmbda,is] = sort(lmbda,'descend') ;
W = W(:,is(1:R)) ;
%clf ; plot(log10(abs(lmbda)),'o')

% Invariance
Wup = W(2:end-1,:) ;
Wdwn = .5*(W(1:end-2,:)+W(3:end,:)) ;
F = Wup\Wdwn ;
z = eig(F) ;

% Wavenumber estimation
k = acos(z) ;

% Amplitude estimation
% uses A.cos(kx+phi) = Up.exp(ikx)+Um.exp(-ikx)
% with Up = A/2.exp(iphi) and Um = A/2.exp(-iphi)
% model matrix
x = (0:nX-1)' ;
V = exp(1i*[k;-k].'.*x) ;
% conditionning
v0 = max(abs(V),[],1) ;
V = V*diag(1./v0) ;
% estimation
U = V\S ;
U = diag(1./v0)*U ; 
% amplitudes & phases
Up = U(1:R,:) ; Um = U(R+1:end,:) ;
PHI = -1i*log(Up./Um)/2 ; 
A = Up.*exp(-1i*PHI) + Um.*exp(1i*PHI) ;


