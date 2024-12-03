function [k,U] = lsfit(S,k,fun)
% Least-square fitting of a signal S with a sum of exponentials/cosines
% Signal model:
%   (exp): s(x,y) = sum_r^R{ a_r(y) exp(1i*k_r*x) }
%   (cos): s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
if nargin<3 ; fun = 'exp' ; end
[nX,nY,nF] = size(S,1:3) ;
[R,~] = size(k) ;

tolK = 1e-5 ; sqrt(eps) ;
maxIt = 300 ;
alphanorm = 2 ;
verbose = false ;
step = 1 ;


% Initialize
    x = (0:nX-1)' ;
    switch fun
        case 'exp'
            U = NaN([R nY nF]) ; % amplitudes
            %Sc = NaN([nX nY nF R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
            Jk = eye(R) ; % wavenumber selection
        case 'cos'
            U = NaN([2*R nY nF]) ; % amplitudes
            %Sc = NaN([nX nY nF 2*R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
            Jk = [eye(R);-eye(R)] ; % wavenumber selection
    end
    if R==0 ; return ; end
    
% Optimization loop
    w = 1 ;
    for ff = 1:nF
        du = inf ; dk = inf ; it = 0 ;
        while it<maxIt && any(abs(dk)>tolK,'all')
            it = it+1 ;
        % Estimate amplitudes U
            % model matrix
                V = exp(1i*(Jk*k(:,ff)).'.*x) ;
            % conditionning
                v0 = max(abs(V),[],1) ;
                Vc = V*diag(1./v0) ;
            % estimation
                u = (Vc'*Vc)\(Vc'*S(:,:,ff)) ;
                u = diag(1./v0)*u ; % inverse conditionning
            % update
                du = u-U(:,:,ff) ;
                U(:,:,ff) = u ; 
        % Update poles k
            sc = permute(V,[1 3 2]).*permute(U(:,:,ff),[3 2 1]) ; % signal components [nX nY (1or2)*R]
            sm = sum(reshape(sc,[],size(sc,3)),2) ; % signal model [nX*nY 1] ;
            ds_dk = 1i*reshape(x.*sc,[],size(sc,3))*Jk ; % signal derivatives [nX*nY R]
            r = sm-reshape(S(:,:,ff),[],1) ; % residue
            w = abs(r).^(alphanorm-2) ;
            j = ds_dk'*(w.*r) ; % jacobian
            H = ds_dk'*(w.*ds_dk) ; % Hessian
            dk = -H\j ;  % update
            k(:,ff) = k(:,ff) + step*dk ;
        % Display infos
            if verbose
                disp("LSFIT" ...
                        + " | it: " + string(it) ...
                        + " | max(abs(dk)) =  " + string(max(abs(dk(:)))) ...
                    ) ;
            end
        end
    end
end