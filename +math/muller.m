function xip1 = muller(fun,x0,maxIt,verbose)
% MULLER Find the zeros of a function via the muller algorithm
tolX = max(sqrt(eps),1e-6*abs(x0)) ;
tolF = sqrt(eps) ; 
if nargin<3 || isempty(maxIt) ; maxIt = 100 ; end
if nargin<4 || isempty(verbose) ; verbose = false ; end

% Initialization with tree points
if iscell(x0) % the tree first points are given
    xip1 = x0{end} ; xi = x0{end-1} ; xim1 = x0{end-2} ;
else % build points around (not aligned !)
    deltaX = 1e1*tolX ;
    xip1 = x0 ; xi = x0+deltaX ; xim1 = x0-3*(1+1i)*deltaX ;
end
fi = fun(xi) ; fim1 = fun(xim1) ;

%tolF = max(tolF,1e-9*fi) ;

% another test: see https://kilyos.ee.bilkent.edu.tr/~microwave/programs/utilities/numeric1/infoMuller.htm
it = 0 ; 
converged = false(size(x0)) ;
while it<maxIt
% Swap the points & functions
    [xi,xim1,xim2] = deal(xip1,xi,xim1) ;
    [fi,fim1,fim2] = deal(fun(xi),fi,fim1) ;
% Updating scheme
    q = (xi-xim1)./(xim1-xim2) ;
    A = q.*fi - q.*(1+q).*fim1 + q.^2.*fim2 ;
    B = (2*q+1).*fi - (1+q).^2.*fim1 + q.^2.*fim2 ;
    C = (1+q).*fi ;
    sqrtB2m4AC = sqrt(B.^2-4*A.*C) ;
    % Find the closest solution to x0
    Dp = B+sqrtB2m4AC ; Dm = B-sqrtB2m4AC ;
    DpgtDm = abs(Dp)>abs(Dm) ;
    D = Dp.*DpgtDm + Dm.*(~DpgtDm) ;
    % Prevent the zero-curvature case & switch to secant method
    isAsmall = abs(A./fi)<tolF ;
    D(isAsmall) = B(isAsmall) ;
    dx = - 2*(xi-xim1).*(C./D) ;
% Prevent NaNs..
    dx(isnan(dx)) = tolX(isnan(dx)) ;
% New guess point
    xip1(~converged) = xi(~converged) + dx(~converged) ;
% Convergence of the point ?
    converged = converged | abs(dx)<tolX ; 
% Convergence of the function ?
    converged = converged | abs(fi)<tolF ; 
% iteration n°
    it = it+1 ;
    if verbose
        disp("MULLER"...
                + " | it: "+string(it)...
                + " | converged: "+string(sum(converged(:)))+"/" + string(numel(x0)) ...
                ...+ " | max(dx) = "+string(max(abs(dx(~converged)))) ...
            ) ;
    end
% Break the loop ?
    if all(converged(:)) ; break ; end
end

end

