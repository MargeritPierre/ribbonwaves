classdef poly
%POLY implements polynom manipulation methods
properties
    coeffs % [nCoeffs nPoly]
end

%% BASE METHODS & PROPERTIES
methods
    function this = poly(c,n)
    % POLY Construct the polynomial representation
        if isa(c,'math.poly') ; this = c ; return ; end
        %if isvector(c) ; c = c(:) ; end
        if nargin<2 ; n = (0:size(c,1)-1)' ; end
        %if isvector(n) ; n = n(:) ; end
        % Build the polynomial 
        n = n + zeros(size(c)) ; % orders
        p = repmat(1:size(c,2),[size(c,1) 1]) ; % polynom index
        this.coeffs = full(sparse(n+1,p,c)) ;
    end
    
    function y = eval(P,x)
    % Evaluate the polynomial at given coordinates x
        X = (x(:).^(0:P.N-1)) ; 
        y = reshape(X*P.coeffs,[size(x) size(P.coeffs,2)]) ;
    end
    
    function N = N(P) ; N = size(P.coeffs,1) ; end
    function nP = npoly(P) ; nP = size(P.coeffs,2) ; end
end


%% OPERATIONS ON POLYNOMIALS
methods
    function P = nnz(P)
    % Keep only non-zero leading coefficients
        lead = find(any(P.coeffs~=0,2),1,'last') ; % leading order
        if isempty(lead) ; lead = 1 ; end
        P.coeffs = P.coeffs(1:lead,:) ;
    end
    
    function P = setorder(P,N)
    % Truncate a polynomial at a given order
        if N<P.N
            P.coeffs = P.coeffs(1:N,:) ;
        else
            P.coeffs(P.N+1:N,:) = 0 ;
        end
    end
    
    function P = keeporders(P,q)
    % Build a new polynomial so that P(y) = P(x^q)
        P.coeffs = P.coeffs(1:q:end,:) ;
    end

    function Pc = plus(Pa,Pb)
    % Sum of two polynomials
        [Pa,Pb] = deal(math.poly(Pa),math.poly(Pb)) ;
        n = max(Pa.N,Pb.N) ;
        Pa.coeffs(end+1:n,:) = 0 ; Pb.coeffs(end+1:n,:) = 0 ; 
        Pc = math.poly(Pa.coeffs+Pb.coeffs) ;
    end

    function Pc = minus(Pa,Pb)
    % Substraction of two polynomials
        Pc = math.poly(Pa) + (-math.poly(Pb)) ;
    end
    
    function P = uplus(P) ; end % P = + P : neutral
    function P = uminus(P) ; P.coeffs = - P.coeffs ; end % P = - P : opposite
    
    function Pc = mtimes(Pa,Pb,N)
    % Product of two polynomials
        [Pa,Pb] = deal(math.poly(Pa),math.poly(Pb)) ;
        a = Pa.coeffs ; b = Pb.coeffs ;
        if Pa.npoly==1 || Pb.npoly==1
            c = convn(a,b) ;
        else
            Nc = Pa.N + Pb.N - 1 ;
            nP = Pa.npoly ; % ==Pb.npoly
            c = zeros(Nc,nP) ;
            for pp = 1:nP
                c(:,pp) = conv(a(:,pp),b(:,pp)) ; 
            end
        end
        Pc = math.poly(c) ;
        if nargin>=3 && ~isempty(N) ; Pc.coeffs = Pc.coeffs(1:N,:) ; end
    end
    
    function Pc = times(varargin) 
    % Same as mtimes
        Pc = mtimes(varargin{1},varargin{2:end}) ;
    end
    
    function [Q,R] = mrdivide(A,B)
    % Polynomial division
        [A,B] = deal(math.poly(A),math.poly(B)) ;
        [a,b] = deal(A.coeffs,B.coeffs) ;
        nPoly = max(size(a,2),size(b,2)) ;
        [a,b] = deal(a+zeros(1,nPoly),b+zeros(1,nPoly)) ;
        q = zeros(A.N-B.N+1,nPoly) ;
        r = zeros(A.N,nPoly) ;
        for pp = 1:nPoly 
            [q(:,pp),r(:,pp)] = deconv(a(:,pp),b(:,pp)) ;
        end
        [Q,R] = deal(math.poly(q),math.poly(r)) ;
    end
    
    function varargout = rdivide(varargin)
    % Same as mrdivide
        varargout = cell(1,nargout) ;
        [varargout{:}] = mrdivide(varargin{:}) ;
    end
    
    function Pc = comp(Pa,Pb,N)
    % Polynomial composition so that C(x) = A(B(x))
    % uses Horner's algorithm. see https://ens-lyon.hal.science/ensl-00546102/document
        [Pa,Pb] = deal(math.poly(Pa),math.poly(Pb)) ;
        [Pa,Pb] = deal(Pa.nnz,Pb.nnz) ; % avoids exponential growth of the polynoms
        Pc = math.poly(0) ; % start with zero polynomial
        for nn = Pa.N:-1:1
            Pc = Pc*Pb + Pa.coeffs(nn,:) ;
            if nargin>2 ; Pc = Pc.setorder(min(Pc.N,N)) ; end
        end
    end
    
    function P = power(P,a) ; P = P.mpower(a) ; end
    function P = mpower(P,a) 
    % Q(x) = P(x)^a
        if mod(a,1)==0 % integer powers
            P = math.poly(full(sparse(a+1,1,1))).comp(P) ;
        else % non-integer powers
            P = P.pow(a) ; 
        end
    end

    function Pd = diff(P,d)
    % d-th derivative of a polynomial
        P = math.poly(P) ;
        n = (d:P.N-1)' ;
        Pd = math.poly(factorial(n)./factorial(n-d).*P.coeffs(d+1:end,:)) ;
    end
    
    function P = ctranspose(P)
    % Prime (') as the derivation sign =)
        P = diff(P,1) ;
    end
    
    function r = roots(P)
    % compute the roots of the polynomial
        if size(P.coeffs,1)==1 
            r = roots(flip(P.coeffs)) ;
        else
            c = num2cell(P.coeffs,1) ;
            r = cellfun(@(c)roots(flip(c)),c,'uni',false) ;
            N = cellfun(@numel,r) ;
            if min(N)~=max(N) % not the same number of roots evrywhere
                N = P.N-1 ; 
                r = cellfun(@(r)[r;NaN(N-numel(r),1)],r,'uni',false) ;
            end
            r = cat(2,r{:}) ;
        end
    end
end

%% ASSIGNMENT BEHAVIOR
methods
    function varargout = subsref(P,s)
    % Referencing the object
        switch s(1).type % P(...)
            case '()'
                x = s.subs{1} ;
                if isa(x,'math.poly') % P(Q) does the polynomial composition
                    varargout = {P.comp(x)} ;
                else % P(x) evaluates the polynomial
                    varargout = {P.eval(x)} ;
                end
            case '{}' % select a set of coefficients/polynoms in the list
                s(1).type = '()' ; 
            % P{i} extracts the full i-th polynom(s) by default
                if numel(s(1).subs)==1 ; s(1).subs = {':' s(1).subs{1}} ; end
            % extract the coefficients
                P = math.poly(builtin('subsref',P.coeffs,s(1))) ;
            % If other subref mechanism is involved..
                if numel(s)>1 ; P = P.subsref(s(2:end)) ; end
                varargout = {P} ;
            otherwise % P{...}, P.(..) etc
                [varargout{1:nargout}] = builtin('subsref',P,s) ;
        end
    end
    
    function P = horzcat(varargin)
    % Concatenate polynomials
        N = max(cellfun(@(P)P.N,varargin)) ;
        varargin = cellfun(@(P)P.setorder(N),varargin,'uni',false) ;
        varargin = cellfun(@(P)P.coeffs,varargin,'uni',false) ;
        P = math.poly(cat(2,varargin{:})) ;
    end
end



%% TAYLOR EXPANSIONS OF USUAL FUNCTIONS
% see https://fr.wikipedia.org/wiki/S%C3%A9rie_de_Taylor
methods (Static)
    function P = X(N) ; P = math.poly([0;1;zeros(max(0,N-2),1)]) ; end % identity P(x)= x
    function P = expX(N) ; n = (0:N-1)' ; P = math.poly(1./factorial(n),n) ; end
    function P = cosX(N) ; n = (0:2:N-1)' ; P = math.poly((-1).^(n/2)./factorial(n),n) ; end
    function P = sinX(N) ; n = (1:2:N-1)' ; P = math.poly((-1).^((n-1)/2)./factorial(n),n) ; end
    function P = sincX(N) ; n = (1:2:N)' ; P = math.poly((-1).^((n-1)/2)./factorial(n),n-1) ; end
    function P = coshX(N) ; n = (0:2:N-1)' ; P = math.poly(1./factorial(n),n) ; end
    function P = sinhX(N) ; n = (1:2:N-1)' ; P = math.poly(1./factorial(n),n) ; end
    function P = log1pX(N) ; n = (1:N-1)' ; P = math.poly((-1).^(n+1)./n,n) ; end % log(1+x) 
    function P = logX(N) ; P = comp(math.poly.log1pX(N),[-1;1],N) ; end % log(x) = log(1+(-1+x))
    function P = pow1pX(N,a) % (1+x)^a
        n = (1:N-1)' ;
        if mod(a,1)==0 % integer powers
            a_n = factorial(a)./factorial(n)./factorial(a-n) ; % binomial coeff
        else % non-integer powers
            a_n = gamma(a+1)./gamma(n+1)./gamma(a-n+1) ; % generalized binomial coeff
        end
        P = math.poly([1;a_n]) ; 
    end
    function P = powX(N,a) 
        if mod(a,1)==0 % integer powers
            P = math.poly(sparse(a+1,1,1,N,1)) ;
        else % (x)^a = (1+(-1+x))^a
            P = comp(math.poly.pow1pX(N,a),[-1;1],N) ; 
        end
    end
    function P = inv1mX(N) ; P = math.poly(ones(N,1)) ; end % 1/(1-x) 
    function P = taylor(f,N,x0)
    % Compute the taylor expansion of the function f 
    % at order N-1
    % around the coordinate x0
        if nargin<3 ; x0 = 0 ; end
    % Use symbolic math to compute the derivatives
        syms x
        an = cell(N,1) ;
        df = f(x) ; % zero derivative order
        for n = 1:N
            an{n} = double(subs(df,x,x0)) ;
            df = diff(df) ; % next derivative order
        end
        an = cat(1,an{:}) ;
        an = an./(factorial(0:N-1)') ;
        P = math.poly(an) ; % P = sum_{ (f^(n)(x0)/n!) X^n}
        P = P.comp([-x0;ones(size(x0))]) ; % P = sum_{ (f^(n)(x0)/n!) (X-x0)^n}
    end
end

%% FUNCTIONS OF A POLYNOMIAL...
methods
    function Q = exp(P) ; Q = math.poly.expX(P.N).comp(P) ; end
    function Q = cos(P) ; Q = math.poly.cosX(P.N).comp(P) ; end
    function Q = sin(P) ; Q = math.poly.sinX(P.N).comp(P) ; end
    function Q = sinc(P) ; Q = math.poly.sincX(P.N).comp(P) ; end
    function Q = cosh(P) ; Q = math.poly.coshX(P.N).comp(P) ; end
    function Q = sinh(P) ; Q = math.poly.sinhX(P.N).comp(P) ; end
    function Q = log1pP(P) ; Q = math.poly.log1pX(P.N).comp(P) ; end % log(1+P) 
    function Q = log(P) ; Q = math.poly.logX(P.N).comp(P) ; end % log(P) = log(1+(-1+P))
    function Q = pow1pP(P,a) ; Q = math.poly.pow1pX(P.N,a).comp(P) ; end % (1+P)^a
    function Q = pow(P,a) % P^a
        if mod(a,1)==0 % integer powers
            Q = math.poly.powX(P.N,a).comp(P) ;
        else % non-integer a: (P)^a = sc^a*(1+(-1+P/sc))^a
        % Scaling factor sc
            [~,is] = max(P.coeffs~=0,[],1) ;
            sc = P.coeffs(sub2ind(size(P.coeffs),is,1:P.npoly)) ;
        % Arrange the scaling
            sc = abs(sc) ;
            sc(sc<1) = 1 ;
        % Apply
            Q = (sc.^a).*math.poly.pow1pX(P.N,a).comp((P./sc)-1) ;  
        end
    end
    function Q = sqrt(P) ; Q = pow(P,.5) ; end
    function Q = inv1mP(P) ; Q = math.poly.inv1mX(P.N).comp(P) ; end
end


%% DISPLAY BEHAVIOR
methods
    function disp(P)
    % Control how a polynomial will be displayed
        builtin('disp',P)
    % Add the polynomial expansion to the display
        P = P.nnz ;
        idx = (1:P.N-1)' ;
        Xstr = ["";strcat(" X^",string(idx))] ;
        Cstr = cellfun(@mat2str,num2cell(P.coeffs,2),'uni',false) ;
        str = strcat(Cstr,Xstr) ;
        delimiter = " + " ;
        if P.npoly>1 ; delimiter = newline + delimiter ; end
        str = strjoin(str,delimiter) ;
        disp(str)
%         fprintf('%s',str) ;
%         fprintf('%s',newline) ;
    end
end











%% UNITARY TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
    function tests
    %% Usual operations
        clc
        P = math.poly([-1;1])  % construction
        P = +P*-P + P*[0;0;2] - P'' % operations
        P(math.poly([0;1]))-P % composition with x is the same polynom
        r = P.roots
        Pr = P(r) % evaluate the polynom at its roots
        
        [Q,R] = P/[1 2;-1 1] %polynomial division and residue
        
        
    %% Function approximation
        x = 3*linspace(-1,1,1000)' ;
        N = 200 ;
        X = math.poly.X(N) ;
        switch 11
        % Taylor series (static methods)
            case 1 ; f = exp(x) ; P = math.poly.expX(N) ;
            case 2 ; f = cos(x) ; P = math.poly.cosX(N) ;
            case 3 ; f = sin(x) ; P = math.poly.sinX(N) ;
            case 4 ; f = sin(x)./x ; P = math.poly.sincX(N) ;
            case 5 ; f = cosh(x) ; P = math.poly.coshX(N) ; 
            case 6 ; f = sinh(x) ; P = math.poly.sinhX(N) ;
            case 7 ; f = log(1+x) ; P = math.poly.log1pX(N) ; 
            case 8 ; f = log(x) ; P = math.poly.logX(N) ;
            case 9 ; a = 3/10 ; f = (1+x).^a ; P = math.poly.pow1pX(N,a) ;
            case 10 ; a = 15/10 ; f = (x).^a ; P = math.poly.powX(N,a) ;
            case 11 ; f = 1./(1-x) ; P = math.poly.inv1mX(N) ;
            case 19 ; a = .5 ; y = 10 ; f = (x + y).^a ; P = (X+y).^a ;
        % Function composition
            case 12 ; f = sqrt(1-sin(x).^2) ; P = {sqrt(1-sin(X)^2) cos(X)}  ;
            case 13 ; f = x ; P = log(exp(X)) ; % P = exp(log(Px)) has a very small convergence radius !!!
            case 14 ; a = 9/10 ; f = x ; P = pow(pow(X,a),1/a) ; % not working well
            case 15 ; a = 1/1.5 ; f = x ; P = pow(pow1pP(X,a),1/a)-1 ; % works well !
        % Comparison
            case 16 ; a = 1/pi ; f = x.^a ; P = {pow1pP(X-1,a) , X^a} ;
            case 17 ; a = 1/pi ; f = (1+x).^a ; P = {pow1pP(X,a) (1+X)^a} ;
        % Automatic taylor expansion
            case 18 ; fun = @(x)sin(x) - cos([2 5].*x) + x.^2 + [1 2] ; f = fun(x(:)) ; P = math.poly.taylor(fun,N) ;
        end
        clf ; 
        plot(x,f) ; 
        drawnow ; set(gca,'ylimmode','manual') ; 
        if iscell(P)
            cellfun(@(P)plot(x,P(x),':'),P,'uni',false) ;
        else
            plot(x,reshape(P(x),[],P.npoly),':')
        end
        
    %% Multi-Polynomial 
        N = 23 ;
        theta = math.poly.X(N) ;
        a = 2 ; b = 1 ; phi = pi/5 ;
        X = [a*cos(theta) b*sin(theta)] ;
        X = [X{1}*cos(phi)-X{2}*sin(phi) X{1}*sin(phi)+X{2}*cos(phi)] ;
        x = X(linspace(0,2*pi,100)') ; 
        clf ; axis equal ;  plot(x(:,1),x(:,2)) ;
        
    end
end


end
