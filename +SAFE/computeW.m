function [U,w] = computeW(MESH,MAT,k,nModes)
% COMPUTE THE FREQUENCIES ASSOCIATED TO WAVES IN A 1D WAVEGUIDE 
solver = 'eigs' ; % 'eigs' or ' or 'lobpcg'

% Inverse iteration shift
hmax = norm(range(MESH.X,1)) ;
fc = .5*sqrt(MAT.G/MAT.rho)/hmax ;
sigma = -0.01*(2*pi*fc)^2 ; % eigenvalue shift
    
% Build the interpolation matrices
    MESH = safe.mesh.quad8(MESH) ; % serendipity element
    MESH = safe.mesh.quadrature(MESH,'quad2') ; % generate a quadrature rule
    MESH = safe.mesh.interpolation(MESH) ; % generate interpolation matrices
                         
% Build the SAFE Matrices
    [K0,K1,K2,M] = safe.matrices(MESH,MAT) ;
  
% Zero-K modes
    [U0,w02] = eigs(K0,M,nModes,sigma) ;
    f0 = sqrt(diag(w02))/2/pi ;

% Wavenumber sweep
w = NaN(numel(k),nModes) ;
U = NaN([size(U0) numel(k)]) ;
for kk = 1:numel(k)
    K = K0 + k(kk)*K1 + k(kk)^2*K2 ;
    switch solver
        case 'eigs'
            [U(:,:,kk),w2] = eigs(K,M,nModes,sigma) ;
            w2 = diag(w2) ;
        otherwise
        % Initial basis
            if kk==1 
                Ui = U0 ;
            else 
                Ui = U(:,:,kk-1) ;
            end
        % Iterate
            switch solver
                case 'lobpcg'
                    [U(:,:,kk),w2] = lobpcg(K,M,Ui,sigma) ;
                case 'subspace'
                    [U(:,:,kk),w2] = subspaceIteration(K,M,Ui,sigma) ;
            end
    end
    w(kk,:) = sqrt(w2(:).') ;
end

end





%% UTILS
function [U,lmbda] = subspaceIteration(A,B,U,sigma,tol,maxIt)
    if nargin<4 ; sigma = 0 ; else ; A = A-sigma*B ; end % Inverse iteration shift
    if nargin<5 ; tol = 1e-4 ; end
    if nargin<6 ; maxIt = 100 ; end
   
    nModes = size(U,2) ;
    nModesIter = min(nModes+8,floor(2*nModes)) ;
    
    U(:,end+1:nModesIter) = randn(size(U,1),nModesIter-nModes) ;
    U = U./sqrt(sum((B*U).*conj(U),1)) ;
    
    lmbdaI = Inf ; it = 0 ;

    while it<maxIt 
        it = it+1 ; 
        
        V = B*U ;
        W = A\V ;
        X = B*W ;

        [Q,lmbda] = eig(W'*V,W'*X,'vector') ;
        [lmbda,is] = sort(lmbda,'ascend') ;
        lmbda = lmbda(1:nModes) ;

        U = W*Q(:,is) ;
        U = U./sqrt(sum((B*U).*conj(U),1)) ;

        if max(abs(lmbdaI-lmbda)./abs(lmbda))<=tol ; break ; end
        lmbdaI = lmbda ;
    end
    lmbda = lmbda + sigma ;
    U = U(:,1:nModes) ;
end

function [U,lmbda] = lobpcg(A,B,U,sigma,tol,maxIt)
    if nargin<4 ; sigma = 0 ; else ; A = A-sigma*B ; end % Inverse iteration shift
    if nargin<5 ; tol = 1e-4 ; end
    if nargin<6 ; maxIt = 100 ; end
    
    nModes = size(U,2) ;
    nModesIter = nModes ; min(nModes+8,floor(2*nModes)) ;
    
    U(:,end+1:nModesIter) = randn(size(U,1),nModesIter-nModes) ;
    
    U = U./sqrt(abs(sum((B*U).*conj(U),1))) ;
    lmbdaI = sum((A*U).*conj(U),1) ; 
    
    it = 0 ; Um = U ;
    P = zeros(size(U,1),0) ;

    while it<maxIt 
        it = it+1 ; 
        
        %R = A*U - (B*U).*lmbdaI(:).' ; % residual
        %R = A\(A*U - (B*U).*lmbdaI(:).') ; % precond. residual
        R = U - A\(B*U).*lmbdaI(:).' ; % precond. residual
        
%         R = R - (U*((B*U)'*R)) ; % R perpto U
%         RBR = R'*B*R ; 
%         RBR = .5*(RBR+RBR') ; 
%         RBR = chol(RBR) ;
%         R = R/RBR ; % R orthonormal
        
        if it>1 
            P = Um-U*(U\Um) ;
            Um = U ;
        end
        X = [U R P] ;
        
        Ax = X'*(A*X) ; 
        Bx = X'*(B*X) ;
        
        [Q,rho] = eig(Ax,Bx,'vector','chol') ;
        [~,is] = sort(rho,'ascend') ;
        U = X*Q(:,is(1:nModesIter)) ;
        
        U = U./sqrt(abs(sum((B*U).*conj(U),1))) ;
        lmbda = sum((A*U).*conj(U),1) ; 
        
        [lmbda,is] = sort(lmbda,'ascend','comparisonmethod','real') ;
        U = U(:,is) ;

        %clf ; plot(real(lmbda),'o') ; drawnow
        if max(abs(lmbdaI(1:nModes)-lmbda(1:nModes))./abs(lmbda(1:nModes)))<=tol 
            break ; 
        end
        lmbdaI = lmbda ;
    end
    lmbda = lmbda(1:nModes) + sigma ;
    U = U(:,1:nModes) ;
end