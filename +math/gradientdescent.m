function x = gradientdescent(f,x,maxIt)
%GRADIENTDESCEND Find the zero of a function via gradient descent
tolX = max(sqrt(eps),1e-6*abs(x)) ;
tolF = sqrt(eps) ; 
if nargin<3 ; maxIt = 100 ; end
deltaX = 1e1*tolX ;
step = .5 ;

it = 0 ;
while it<maxIt
    df_dx = .5*(f(x+deltaX)-f(x-deltaX))./deltaX ;
    fx = f(x) ;
    dx = - (conj(df_dx(:)).*fx(:))./(abs(df_dx(:)).^2) ;
    x(:) = x(:) + step*dx ;
    it = it+1 ; 
end

end

