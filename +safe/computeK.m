function [U,k] = computeK(MESH,MAT,w,nModes,redux)
% COMPUTE THE FREQUENCIES ASSOCIATED TO WAVES IN A 1D WAVEGUIDE 
if nargin<5 ; redux = false ; end

% Inverse iteration shift
    sigma = 0 ;
                         
% Build the SAFE Matrices
    [K0,K1,K2,M] = safe.matrices(MESH,MAT) ;
    O = sparse(size(M,1),size(M,2)) ;
    I = speye(size(M)) ;
    
% Parametric eigenvalue problem
    A = @(p)blkdiag(p^2*M-K0,I) ;
    B = [K1 K2 ; I O] ;
    [U,k] = safe.parameig(A,B,w,nModes,sigma,redux) ;
    U = U(1:end/2,:,:) ;

% Wavenumber sweep
%     k = NaN(nModes,numel(w)) ;
%     U = NaN([size(M,1) nModes numel(w)]) ;
%     wtbr = waitbar(0,'Compute K...') ; st = tic ;
%     for ww = 1:numel(w)
%         A = blkdiag(w(ww)^2*M-K0,I) ; 
%         B = [K1 K2 ; I O] ;
%         [UkU,K] = eigs(A,B,nModes,sigma) ;
%         k(:,ww) = diag(K) ;
%         U(:,:,ww) = UkU(1:end/2,:) ;
%         if toc(st)>.25 ; wtbr = waitbar(ww/numel(w),wtbr) ; st = tic ; end
%     end
%     delete(wtbr) ;

end
