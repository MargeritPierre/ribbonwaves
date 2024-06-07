function vr = rayleigh(nu,vt,approx)
% RAYLEIGH computes the rayleigh wave velocity
if nargin<2 || isempty(vt) ; vt = 1 ; end
if nargin<3 ; approx = numel(nu)>1 ; end

if approx % use pre-computed polynomial approximation (see below)
    a = [0.874031236600605;0.195420648683553;-0.0402212021593556;-0.0681876727754409;0.0121293927301497;0.0500284821714939;0.00598542571353657;-0.0272004654530807;-0.0129775670388534] ;
    p = 0:numel(a)-1 ;
    A = nu(:).^p ; 
    vr = reshape(A*a,size(nu)) ;
else % use exact solution
% Initial guess with approximation
    vr0 = analytical.rayleigh(nu,1,true) ;
    beta = vr0.^2 ; 
    alpha = .5*(1-2*nu)./(1-nu) ;
% Find the root of the characteristic equation
    for bb = 1:numel(beta)
        Fb = @(beta) (beta-2).^4 - 16*(beta-1).*(alpha(bb).*beta-1) ;
        beta(bb) = fzero(Fb,beta(bb)) ;
    end
% Normalized Velocity
    vr = sqrt(beta) ;
end

% Final velocity
    vr = vr.*vt ;
end




function coeffs
%% ESTIMATE THE VELOCITY APPROXIMATION POLYNOMIAL


ord = 8 ; % polynomial order
nu = linspace(-1+eps,.5-eps,1000) ;
vr = analytical.rayleigh(nu,1,false) ;
vra = analytical.rayleigh(nu,1,true) ;

p = 0:ord ;
A = nu(:).^p ; 
a = A\vr(:) ;
vrp = A*a ;

clf ; 
plot(nu,vr) ;
plot(nu,vra,'-.') ;
plot(nu,vrp,':') ;

disp("a = " + mat2str(a) + " ;")
err = norm(vrp-vr(:))/norm(vr)


end