function [k,A,PHI] = esprit(W,fun,S)
% PARAMETRIC ESTIMATION OF A COMBINATION OF EXP|COS FUNCTIONS
% Signal model :
%   (exp): s(x,y) = sum_r^R{ a_r(y) exp(1i*k_r*x) }
%   (cos): s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
if nargin<2 ; fun = 'exp' ; end
[m,R,nF] = size(W,1:3) ;

% Shifted subspaces
switch fun
    case 'exp'
        Wup = W(1:end-1,:,:) ;
        Wdwn = W(2:end,:,:) ;
%         Wup = permute(W(1:end-1,:,:),[2 4 3 1]) ; % [R 1 nF nX]
%         Wdwn = permute(W(2:end,:,:),[4 2 3 1]) ; % [1 R nF nX]
        % we = permute(W(end,:,:),[2 1 3]) ;
        % iw2 = 1./(1-sum(abs(we).^2,2)) ;
        % WW = sum(conj(Wup).*Wdwn,4) ;
        % F = WW + (iw2.*we).*sum(we.*WW,1) ;
    case 'cos'
        Wup = W(2:end-1,:,:) ;
        Wdwn = .5*(W(1:end-2,:,:)+W(3:end,:,:)) ;
end

% Spectral matrices
F = NaN(R,R,nF) ;
for ff = 1:nF
    F(:,:,ff) = Wup(:,:,ff)\Wdwn(:,:,ff) ;
end

% Pole matrices
z = NaN(R,nF) ;
for ff = 1:nF
    z(:,ff) = eig(F(:,:,ff)) ;
end

% wavenumbers
switch fun
    case 'exp'
        k = -1i*log(z) ;
    case 'cos'
        k = acos(z) ;
end




return
%%

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


