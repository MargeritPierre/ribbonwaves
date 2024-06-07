function [U,lmbda] = parameig(A,B,p,n,sigma,redux)
% Solve a paramteric eigenvalue problem [A(p)-lmbda.B(p)].U=0
% A and B are square (opt. sparse) matrices 
%   or function handles of the parameter p
if nargin<6 ; redux = false ; end

% Reduction parameters
    pr = p(round(linspace(1,end,5))) ; % number of 'key' parameters
    pr(end+1) = 0 ; p(abs(p)==min(abs(p))) ;
    r = 2*n ; min(n+8,2*n) ; % extended basis for reduction
    tolR = sqrt(eps) ; % tolerance to cull basis vectors

% Size info
    if isa(A,'function_handle') ; Ap = A(p(1)) ; else ; Ap = A ; end
    if isa(B,'function_handle') ; Bp = B(p(1)) ; else ; Bp = B ; end
    nDOF = size(Bp,1) ; 
    nP = numel(p) ;

% Model order reduction via coarse parameter sweep
% Compute the reduced basis
    if redux
        nPr = numel(pr) ;
        Vr = zeros(nDOF,r,nPr) ;
        Ar = sparse(0) ; Br = sparse(0) ;
    % Compute the eigenvectors for the key parameter
        wtbr = waitbar(0,'Reduced basis...') ;
        for pp = 1:nPr
        % Evaluate matrix operators if needed
            if isa(A,'function_handle') ; Ap = A(p(pp)) ; end
            if isa(B,'function_handle') ; Bp = B(p(pp)) ; end
        % Compute eigenvectors
            [Vr(:,:,pp),~] = eigs(Ap,Bp,r,sigma) ;
            Ar = Ar + Ap ; % mean A-operator
            Br = Br + Bp ; % mean B-operator
        % Progress bar
            wtbr = waitbar(pp/nPr,wtbr) ;
        end
        delete(wtbr) ;
    % Orthogonalize the basis
        Vr = Vr(:,:) ;
        Ao = (1/nPr)*(Vr'*Ar*Vr) ;
        Bo = (1/nPr)*(Vr'*Bp*Vr) ; 
        [Q,lq] = eig(Ao,Bo,'vector') ;
        overtol = abs(lq)>tolR*max(abs(lq(~isinf(abs(lq))))) ;
        disp("Ratio of modes kept: "+string(sum(overtol)/numel(overtol))) ;
        Q = Q(:,overtol) ;
        Vr = Vr*Q ;
    % Apply the projection
        nDOF = size(Vr,2) ;
        if ~isa(A,'function_handle') ; Ap = Vr'*Ap*Vr ; end
        if ~isa(B,'function_handle') ; Bp = Vr'*Bp*Vr ; end
    end

% Parameter sweep
    lmbda = NaN(n,nP) ;
    U = NaN([nDOF n nP]) ;
    wtbr = waitbar(0,'Parametric EVP...') ; st = tic ;
    for pp = 1:nP
    % Evaluate matrix operators if needed
        if isa(A,'function_handle') ; Ap = A(p(pp)) ; end
        if isa(B,'function_handle') ; Bp = B(p(pp)) ; end
    % Apply projection if needed
        if redux && isa(A,'function_handle') ; Ap = Vr'*Ap*Vr ; end 
        if redux && isa(B,'function_handle') ; Bp = Vr'*Bp*Vr ; end 
    % Compute eigenvalues
        if issparse(Ap) || issparse(Bp)
            [uu,lm] = eigs(Ap,Bp,n,sigma) ;
            lm = diag(lm) ;
        else
            [uu,lm] = eig(Ap,Bp,'vector') ;
            [lm,is] = sort(lm-sigma,'ascend','comparisonmethod','abs') ;
            lm = lm(1:n) ;
            uu = uu(:,is(1:n)) ;
        end
    % Backup
        lmbda(:,pp) = lm ;
        U(:,:,pp) = uu ;
    % Progress bar
        if toc(st)>.25 ; wtbr = waitbar(pp/nP,wtbr) ; st = tic ; end
    end
    delete(wtbr) ;
    
% Back to full basis
    if redux
        U = reshape(Vr*U(:,:),[size(Vr,1) size(U,2:ndims(U))]) ;
    end


end