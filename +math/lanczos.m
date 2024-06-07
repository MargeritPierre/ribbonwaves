function [a,b,V,r] = lanczos(opA,m,r,a,b,V)
% Lanczos partial tridiagonalization algorithm of a symmetric matrix A[n n sz]
% A is given by the anonymous operator opA = @(X)A*X, with X[n p sz]
% sz is a size for vectorization
% m is the queried number of Lanczos vectors
% r is the current residual vector with V'*r = 0 and norm(r) = 0
% V[n k sz] is the (partial) starting basis of vectors
% a[m 1 sz] contains the diagonal elements
% b[m-1 1 sz] contains the first upper/lower (sym) elements
[n,~,sz] = size(r(:,:,:)) ;
if nargin<4 ; a = NaN(0,1,sz) ; b = NaN(0,1,sz) ; V = NaN(n,0,sz) ; end
k = size(V,2) ;
% Memory allocation or restart
    V(:,k+1:m,:) = 0 ;
    a(k+1:m,:) = 0 ;
    b(k+1:m,:) = 0 ;
% Iterations
    for ii = k+1:m
        V(:,ii,:) = r ; % backup
        r = opA(r) ; % new residual
    % Three-term recurrence
        a(ii,:,:) = real(pagemtimes(V(:,ii,:),'ctranspose',r,'none')) ; %=real(sum(V(:,ii,:).*conj(r),1)) ;
        r = r - a(ii,:,:).*V(:,ii,:) ;
        if ii>k+1 ; r = r - b(ii-1,:,:).*V(:,ii-1,:) ; end
    % Re-orthogonalization
        r = r - pagemtimes(V,pagemtimes(V,'ctranspose',r,'none')) ; %= r-V(:,1:ii)*V(:,1:ii)'*r
        b(ii,:,:) = sqrt(sum(abs(r).^2,1)) ;
        r = r./b(ii,:,:) ; % normalized residual
    end
end