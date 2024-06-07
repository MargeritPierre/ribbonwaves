function [edg,id] = adjacency(DIST,tol,unique)
% Build the adjacency matrix associated to "close" points (with a given distance
% measure)
% DIST[N1 N2 szD(3:end)]
if nargin<2 || isempty(tol) ; tol = Inf ; end
if nargin<3 ; unique = false ; end

% If unique adjacency is needed (one-to-one matching)
if unique
    nUnique = min(size(DIST,1:2)) ;
    for mm = 1:size(DIST(:,:,:),3)
        d = DIST(:,:,mm) ;
        % Sort by distance
        [~,is] = sort(d(:)) ;
        % Corresponding "next" and "current" indices
        [nn,cc] = ind2sub(size(d),is) ;
        % browse the list
        ii = 0 ;
        while sum(~isnan(d(:)))>nUnique
            ii = ii+1 ;
            if ~isnan(d(nn(ii),cc(ii)))
                d(nn(ii),:) = NaN ;
                d(:,cc(ii)) = NaN ;
                d(nn(ii),cc(ii)) = DIST(nn(ii),cc(ii),mm) ;
            end
        end
        DIST(:,:,mm) = d ;
    end
end

% Truncate with tolerance
DIST(DIST>tol) = NaN ;

% Build connectivity graph edges
id = find(~isnan(DIST(:))) ;
szD = size(DIST,1:3) ;
[e2,e1,nn] = ind2sub(szD,id) ;
edg = [e1 e2] + [nn-1 nn]*max(szD(1),szD(2)) ;

end






