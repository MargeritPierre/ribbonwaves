function [U,w] = computeW(MESH,MAT,k,nModes,redux)
% COMPUTE THE FREQUENCIES ASSOCIATED TO WAVES IN A 1D WAVEGUIDE 
if nargin<5 ; redux = false ; end

% Inverse iteration shift
    hmax = norm(range(MESH.X,1)) ;
    fc = .5*sqrt(real(MAT.G)/MAT.rho)/hmax ;
    sigma = 0 ; %-0.01*(2*pi*fc)^2 ; % eigenvalue shift
                         
% Build the SAFE Matrices
    [K0,K1,K2,M] = safe.matrices(MESH,MAT) ;
    
% Parametric eigenvalue problem
    A = @(p)K0 + p*K1 + p^2*K2 ;
    B = M ;
    [U,w2] = safe.parameig(A,B,k,nModes,sigma,redux) ;
    w = sqrt(w2) ;

% % Wavenumber sweep
% w = NaN(nModes,numel(k)) ;
% U = NaN([size(M,1) nModes numel(k)]) ;
% wtbr = waitbar(0,'Compute W...') ; st = tic ;
% for kk = 1:numel(k)
%     K = K0 + k(kk)*K1 + k(kk)^2*K2 ;
%     [U(:,:,kk),w2] = eigs(K,M,nModes,sigma) ;
%     w(:,kk) = sqrt(diag(w2)) ;
%     if toc(st)>.25 ; wtbr = waitbar(kk/numel(k),wtbr) ; st = tic ; end
% end
% delete(wtbr) ;

end
