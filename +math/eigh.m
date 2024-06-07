function [Vc,dc] = eigh(opA,n_k_p,tol)
% A vectorized version of the thick-restarted Laczos algorithm
% opA is the operator describing the matrix-vector product
% n is the size of the operator
% k is the number of eigenvalues queried
% p is the number of problems

if nargin<3 ; tol = sqrt(eps) ; end
maxIt = 300 ;
eigMethod = 'eig' ; % 'eig' or 'tridiag-QL'

% Size information
n_k_p(end+1:3) = 1 ;
[n,k,p] = deal(n_k_p(1),n_k_p(2),prod(n_k_p(3:end))) ;
m = max(min(n,20),round(2*k)) ; % number of Lanczos vectors

% Does the operator allows index selection ?
iterate = true(p,1) ; p_i = p ; % keep trace of converged indices
if nargin(opA)>1 
    opAi = @(U)opA(U,iterate) ;
else
    opAi = opA ;
end

% Starting lanczos vector
v0 = randn(n,1,p) ;
v0 = v0./sqrt(sum(abs(v0).^2,1)) ;

% Initial lanczos tridiag
[a,b,V,r] = math.lanczos(opAi,m,v0) ;
c = zeros(m,1,p) ;

% Restarting loop
k0 = k ;
it = 0 ;
dc = NaN(k0,1,p) ; Vc = NaN(n,k0,p) ; % converged basis
while 1
    
    % Ritz values of the current tridiagonal matrix
    switch eigMethod
        case 'eig' % Non-vectorized EIG
            Q = NaN(m,m,p_i) ; d = NaN(m,1,p_i) ;
            for pp=1:p_i
                if any(isnan(a(:,:,pp)),'all') || any(isnan(b(:,:,pp)),'all') ; continue ; end
                Tp = diag(a(:,:,pp),0) + diag(b(1:end-1,:,pp),1) + diag(b(1:end-1,:,pp),-1);
                Tp(k+1,1:k) = c(1:k,:,pp)' ; Tp(1:k,k+1) = c(1:k,:,pp) ;
                [Q(:,:,pp),d(:,:,pp)] = eig(Tp,'vector') ;
            end
        case 'tridiag-QL' % via the QL algorithm
         % Convert to tridiagonal form
            if it==0 % first iteration, already tridiagonal
                at = a ; bt = b ; Vt = ones(1,1,p) ;
            else % apply the full lanczos algorithm
                opT = @(U)a.*U ... % Form the almost-tridiagonal matrix operator
                            + [b(1:end-1,:,:).*U(2:end,:,:);zeros([1 size(U,2:3)])] ...
                            + [zeros([1 size(U,2:3)]);b(1:end-1,:,:).*U(1:end-1,:,:)] ...
                            + [zeros([k size(U,2:3)]);sum(c(1:k,:,:).*U(1:k,:,:),1);zeros([m-k-1 size(U,2:3)])] ...
                            + [c(1:k,:,:).*U(k+1,:,:);zeros([m-k size(U,2:3)])] ...
                            ;
                vr = randn(m,1,p_i) ;
                vr = vr./sqrt(sum(abs(vr).^2,1)) ;
                [at,bt,Vt] = math.lanczos(opT,m,vr) ;
            end
        % Extract eigenvalues
            [d,Q] = math.symtrilQL(at,bt,tol,maxIt) ;
            for pp=1:p_i ; Q(:,:,pp) = Vt(:,:,pp)*Q(:,:,pp) ; end
    end

    % Sort eigenvalues
    [d,is] = sort(d,1,'descend') ;
    for pp=1:p_i ; Q(:,:,pp) = Q(:,is(:,pp),pp) ; end

    % Update eigenvectors
    for pp=1:p_i ; V(:,1:k,pp) = V(:,:,pp)*Q(:,1:k,pp) ; end

    % Residual
    c = permute(b(m,:,:).*Q(m,:,:),[2 1 3]) ;

    % Convergence tests
    notConverged = abs(c(1:k,:,:))>tol.*abs(d(1:k,:,:)) ; % [k 1 p_i]
    nConv = sum(~notConverged,1) ; % [1 1 p_i] ;
   
    % Backup converged values
    if nargin(opA)>1 % Does the operator allows index selection ? 
        dc(:,:,iterate) = d(1:k0,:,:) ;
        Vc(:,:,iterate) = V(:,1:k0,:) ;
        iterated = iterate ;
        iterate(iterate) = nConv<k0 ; % [1 1 p_i] ;
    % Constrain the iteration to the needed indices
        opAi = @(U)opA(U,iterate) ;
        ip = iterate(iterated) ; p_i = sum(ip) ;
        a = a(:,:,ip) ; 
        b = b(:,:,ip) ; 
        c = c(:,:,ip) ; 
        d = d(:,:,ip) ;
        V = V(:,:,ip) ;
        r = r(:,:,ip) ;
    else
        dc(:,:,iterate) = d(1:k0,:,iterate) ;
        Vc(:,:,iterate) = V(:,1:k0,iterate) ;
        iterate(iterate) = nConv(iterate)<k0 ; % [1 1 p_i] ;
    end


    % Display infos
    disp("it:"+string(it)+",overtol:"+string(sum(notConverged(:)))+"/"+string(p*k)) ;

    % Algo stop conditions
    if ~any(iterate(:)) ; break ; end
    it = it + 1 ;
    if it>maxIt ; break ; end

    % extend the basis if needed
    %k = k + min(nconv,floor((m-k)/2)) 

    % Thick Restart with Ritz values
    % see https://sdm.lbl.gov/~kewu/ps/trlan-siam.pdf
    % see https://arxiv.org/abs/1902.02064
    a(1:k,:,:) = d(1:k,:,:) ;
    b(1:k,:,:) = 0 ;
    V = V(:,1:k,:) ;
    [a,b,V,r] = math.lanczos(opAi,m,r,a,b,V) ;
    
end

end