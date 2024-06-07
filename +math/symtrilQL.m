function [d,U] = symtrilQL(d,e,tol,maxIt)
% QL eigendecompositon algorithm applied to a symmetric tridiagonal matrix
% d[m 1 sz] and e[m-1 1 sz] are the elements of the diagonal and upper diagonal
% sz is a size for vectorization
% see https://blogs.mathworks.com/cleve/2015/06/29/dubrulle-creates-a-faster-tridiagonal-qr-algorithm/
if nargin<3 ; tol = (eps) ; end % convergence tolerance
if nargin<4 ; maxIt = 10000 ; end % maximum number of QL iterations
    d = d(:,:,:) ; e = e(:,:,:) ;
    [m,~,sz] = size(d,1:3) ;
    if nargout>1 ; U = repmat(eye(m),[1 1 sz]) ; end
    it = 0 ;  jj = 1 ; % last unconverged eigenpair
    while any(jj<m) && it<maxIt
    % Applied shift: see Algo 4.6 of https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
        ss = .5*(d(jj+1,:,:)-d(jj,:,:)) ;
        shift = d(jj,:,:)-e(jj,:,:).^2./(ss+sign(ss).*sqrt(ss.^2+e(jj,:,:).^2)) ;
    % Initialization
        c = ones(m,1,sz) ; s = c ;
        q = d(m,:,:)-shift ;
        u = zeros(1,1,sz) ;
    % Iterations over the vectors
        for ii = m-1:-1:jj
            h = c(ii+1,:,:).*e(ii,:,:) ;
            p = s(ii+1,:,:).*e(ii,:,:) ;
            e(ii+1,:,:) = sqrt(p.^2 + q.^2) ;
            s(ii,:,:) = p./e(ii+1,:,:) ;
            c(ii,:,:) = q./e(ii+1,:,:) ;
            g = d(ii+1,:,:) - u ;
            t = (d(ii,:,:)-g).*s(ii,:,:) + 2.*c(ii,:,:).*h ;
            u = s(ii,:,:).*t ;
            d(ii+1,:,:) = g + u ;
            q = c(ii,:,:).*t - h ;
        end
        d(jj,:,:) = d(jj,:,:) - u ;
        e(jj,:,:) = q ;
    % Eigenvectors
        if nargout>1 
            for ii = m-1:-1:jj
                U(:,ii:ii+1,:) = U(:,ii,:).*[c(ii,:,:) s(ii,:,:)] + U(:,ii+1,:).*[-s(ii,:,:) c(ii,:,:)] ; % = pagemtimes(U(:,ii:ii+1,:),[c s;-s c]) ;
            end
        end
        % U = diag(s(1:end-1,:,:),1)
    % Convergence check
        it = it+1 ;
        unconverged = abs(e(1:m-1,:,:))>tol*abs(d(1:m-1,:,:)+abs(d(2:m,:,:))) ;
        jj = find(any(unconverged,3),1,'first') ;
    end
end