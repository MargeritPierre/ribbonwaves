function [W,lmbda] = subspace(S,m_r,fun)
% Vectorized signal subspace estimation via subspace iteration
% S [l p sz], with:
%   - l the number of "data" samples
%   - p the number of data "realizations"
%   - sz other dimensions for vectorization
% Convergence parameters (hard-coded for now..)
if isscalar(m_r) ; m = [] ; r = m_r ;  
else ; m = m_r(1) ; r = m_r(2) ; end
if nargin<3 ; fun = 'exp' ; end
relTolL = 1e-12 ; sqrt(eps) ; % tolerance on eigenvalues

% Size info
    sz = size(S,1:max(3,ndims(S))) ;
    [q,p,sz] = deal(sz(1),sz(2),sz(3:end)) ; 

% Fourier transform of the signal once for all
    S = permute(S(:,:,:),[1 4 3 2]) ; % [q 1 sz p]
    Fs = fft(S,[],1) ; % [q 1 sz p]

% Sizes of the Hankel/Toeplitz matrices H/T[m n*p] 
% and equivalent operators
    switch fun
        case 'exp'
            if isempty(m) ; m = floor(q/2) ; end
            opCss = @(U,ii)math.expcovtimes(Fs(:,:,ii,:),U) ;
        case 'cos'
            if isempty(m) ; m = floor(q/3) ; end
            opCss = @(U,ii)math.coscovtimes(Fs(:,:,ii,:),U) ;
    end
    
% Use the Lanczos algo
    n_k_p = [m r sz] ;
    [W,lmbda] = math.eigh(opCss,n_k_p,relTolL) ;

end






