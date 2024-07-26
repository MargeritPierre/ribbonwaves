function [k,U,Sc] = esprit(W,fun,S)
% PARAMETRIC ESTIMATION OF A COMBINATION OF EXP|COS FUNCTIONS
% Signal model :
%   (exp): s(x,y) = sum_r^R{ a_r(y) exp(1i*k_r*x) }
%   (cos): s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
if nargin<2 ; fun = 'exp' ; end
[m,R,nF] = size(W,1:3) ;

% Retrieve eigenvalue information
lmbda = sum(abs(W).^2,1) ;
W = W./sqrt(lmbda) ; % renormalize eigenspace

% Shifted subspaces
switch fun
    case 'exp'
        Wup = W(1:end-1,:,:) ;
        Wdwn = W(2:end,:,:) ;
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

% amplitudes
if nargin<3 || nargout<2 ; return ; end

% Amplitude estimation
% uses U.cos(kx+phi) = Up.exp(ikx)+Um.exp(-ikx)
% with Up = U/2.exp(iphi) and Um = U/2.exp(-iphi)
[nX,nY] = size(S,1:2) ;
switch fun
    case 'exp'
        U = NaN([R nY nF]) ; % amplitudes
        Sc = NaN([nX nY nF R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
    case 'cos'
        U = NaN([2*R nY nF]) ; % amplitudes
        Sc = NaN([nX nY nF 2*R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
end
x = (0:nX-1)' ;
for ff = 1:nF
% model matrix
    switch fun
        case 'exp'
            V = exp(1i*k(:,ff).'.*x) ;
        case 'cos'
            V = exp(1i*[k(:,ff);-k(:,ff)].'.*x) ;
    end
% conditionning
    v0 = max(abs(V),[],1) ;
    V = V*diag(1./v0) ;
% estimation
    U(:,:,ff) = V\S(:,:,ff) ;
% reconstruct signal model
    if nargout>2 ; Sc(:,:,ff,:) = permute(V,[1 3 4 2]).*permute(U(:,:,ff),[3 2 4 1]) ; end
% inverse conditionning
    U(:,:,ff) = diag(1./v0)*U(:,:,ff) ; 
end


