function xip1 = muller(fun,x0,maxIt)
% MULLER Find the zeros of a function via the muller algorithm
tolX = max(sqrt(eps),1e-6*abs(x0)) ;
tolF = sqrt(eps) ; 
if nargin<3 ; maxIt = 100 ; end

% Initialization with tree points
deltaX = 1e1*tolX ;
xip1 = x0 ; xi = x0+deltaX ; xim1 = x0-deltaX ;
fi = fun(xi) ; fim1 = fun(xim1) ;

%tolF = max(tolF,1e-9*fi) ;

% another test: see https://kilyos.ee.bilkent.edu.tr/~microwave/programs/utilities/numeric1/infoMuller.htm
it = 0 ; 
converged = false(size(x0)) ;
while it<maxIt
% Convergence of the function ?
    converged = converged | abs(fi)<tolF ; 
    if all(converged) ; break ; end
% Swap the points & functions
    xim2 = xim1 ; xim1 = xi ; xi = xip1 ;
    fim2 = fim1 ; fim1 = fi ; fi = fun(xi) ;
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
    if all(converged) ; break ; end
% iteration n°
    it = it+1 ;
    disp("MULLER"...
            + " | it: "+string(it)...
            + " | converged: "+string(sum(converged(:)))+"/" + string(numel(x0)) ...
            + " | max(dx)= "+string(max(abs(dx(~converged)))) ...
        ) ;
end
end

