function D = distance(X,Y,normalization) 
% Compute the Distances between point sets
% X[nCoord nPx szX(3:end)] and Y[nCoord nPy szY(3:end)]
% D[nPx nPy max(szX(2:end),szY(2:end))] ;
if nargin<2 || isempty(Y) ; Y = X ; end % compute the auto-MAC by default
if nargin<3 || isempty(normalization) ; normalization = 'none' ; end % no normalization by default

% nCoord Complex coordinates are threated as 2*nCoord real coordinates
if ~isreal(X) || ~isreal(Y)
    X = cat(1,real(X),imag(X)) ;
    Y = cat(1,real(Y),imag(Y)) ;
end

% Size information
nDims = max(ndims(X),ndims(Y)) ;
szX = size(X,1:nDims) ;
szY = size(Y,1:nDims) ;
if szX(1)~=szY(1) ; error('the first dimension of X and Y should match') ; end
szOK = szX==1 | szY==1 | szY==szX ;
if any(~szOK(3:end)) ; error('the sizes of X and Y for dimensions>3 should match or be singleton') ; end 

% Distance
X = permute(X,[2,nDims+1,3:nDims,1]) ;
Y = permute(Y,[nDims+1,2,3:nDims,1]) ;
XY = cat(nDims+2,X+0*Y,Y+0*X) ;

% Normalization
switch normalization
    case 'range' % range over all data
        N = max(XY(:))-min(XY(:)) ;
    case 'ranges' % range over individual coordinates
        N = max(XY,[],[1:nDims,nDims+2])-min(XY,[],[1:nDims,nDims+2]) ;
    case 'relative' % relative to coordinate vector magnitudes
        N = sqrt(.5*sum(abs(XY).^2,nDims+(1:2))) ;
    case 'none'
        N = 1 ;
    otherwise
        error('Wrong argument for normalization') ;
end
D = sqrt(sum(abs(diff(XY./N,1,nDims+2)).^2,nDims+1)) ;

end