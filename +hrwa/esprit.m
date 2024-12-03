function [k,U,Sc] = esprit(W,fun,S)
% PARAMETRIC ESTIMATION OF A COMBINATION OF EXP|COS FUNCTIONS
% Signal model :
%   (exp): s(x,y) = sum_r^R{ a_r(y) exp(1i*k_r*x) }
%   (cos): s(x,y) = sum_r^R{ a_r(y) cos(k_r*x + phi_r(y)) }
if nargin<2 ; fun = 'exp' ; end
[m,R,nF] = size(W,1:3) ;

% Retrieve eigenvalue information
lmbda = sum(abs(W).^2,1) ;
W = W./sqrt(lmbda) ; % renormalize eigenspace

% Shifted subspaces
switch fun
    case 'exp'
        Wup = W(1:end-1,:,:) ;
        Wdwn = W(2:end,:,:) ;
    case 'cos'
        Wup = W(2:end-1,:,:) ;
        Wdwn = .5*(W(1:end-2,:,:)+W(3:end,:,:)) ;
end

% Spectral matrices
F = NaN(R,R,nF) ;
for ff = 1:nF
    F(:,:,ff) = Wup(:,:,ff)\Wdwn(:,:,ff) ;
end

% Spectral matrices
F = NaN(R,R,nF) ;
for ff = 1:nF
    F(:,:,ff) = Wup(:,:,ff)\Wdwn(:,:,ff) ;
end

% Pole matrices
z = NaN(R,nF) ;
for ff = 1:nF
    z(:,ff) = eig(F(:,:,ff)) ;
end

% wavenumbers
switch fun
    case 'exp'
        k = -1i*log(z) ;
    case 'cos'
        k = acos(z) ;
end

% amplitudes
if nargin<3 || nargout<2 ; return ; end

% Amplitude estimation
% uses U.cos(kx+phi) = Up.exp(ikx)+Um.exp(-ikx)
% with Up = U/2.exp(iphi) and Um = U/2.exp(-iphi)
[nX,nY] = size(S,1:2) ;
switch fun
    case 'exp'
        U = NaN([R nY nF]) ; % amplitudes
        Sc = NaN([nX nY nF R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
    case 'cos'
        U = NaN([2*R nY nF]) ; % amplitudes
        Sc = NaN([nX nY nF 2*R]) ; % reconstructed signal model components (S=sum(Sm,4)) ;
end
x = (0:nX-1)' ;
for ff = 1:nF
% model matrix
    switch fun
        case 'exp'
            V = exp(1i*k(:,ff).'.*x) ;
        case 'cos'
            V = exp(1i*[k(:,ff);-k(:,ff)].'.*x) ;
    end
% conditionning
    v0 = max(abs(V),[],1) ;
    V = V*diag(1./v0) ;
% estimation
    U(:,:,ff) = V\S(:,:,ff) ;
% reconstruct signal model
    if nargout>2 ; Sc(:,:,ff,:) = permute(V,[1 3 4 2]).*permute(U(:,:,ff),[3 2 4 1]) ; end
% inverse conditionning
    U(:,:,ff) = diag(1./v0)*U(:,:,ff) ; 
end

return ; 


%% WHAT FOLLOWS IS EXPERIMENTAL !!
    
%% STABILIZATION CRITERION
F = cell(R,1) ; % needed for ESTER criterion
ERR = NaN(R,nF) ;
for r = 1:R
    F{r} = NaN(r,r,nF) ;
    for ff = 1:nF
        F{r}(:,:,ff) = Wup(:,1:r,ff)\Wdwn(:,1:r,ff) ;
        ERR(r,ff) = norm(Wup(:,1:r,ff)*F{r}(:,:,ff)-Wdwn(:,1:r,ff),'fro') ;
    end
end

% Pole matrices
z = cell(R,1) ;
for r = 1:R
    z{r} = NaN(r,nF) ;
    for ff = 1:nF
        z{r}(:,ff) = eig(F{r}(:,:,ff)) ;
    end
end

% wavenumbers
switch fun
    case 'exp'
        k = cellfun(@(z)-1i*log(z),z,'uni',false) ;
    case 'cos'
        k = cellfun(@(z)acos(z),z,'uni',false) ;
end

%% STABILIZATION DIAGRAM
ESTER = min(ERR,[],1)./ERR ;
tolEster = 1/1 ;

r = repmat(repelem(1:R,1:R)',[1 nF]) ;
f = repmat(1:nF,[size(r,1) 1]) ;
K = cat(1,k{:}) ;
clf ; plot3(real(K(:)),f(:),r(:),'.','markersize',2) ;

tol = 1/1*(pi/m) ; Nmin = 15 ;
KU = {} ; FU = {} ; KUf = {} ;
KE = {} ; FE = {} ;
for ff = 1:nF
    ku = K(:,ff) ;
    [~,iall] = uniquetol([real(ku) imag(ku)],tol ...
                            ,'datascale',1 ...
                            ,'byrows',true ...
                            ,'outputallindices',true ...
                            ) ;
    ku = cellfun(@(ii)median(ku(ii)),iall) ;
    ku(cellfun(@numel,iall)<Nmin) = [] ;
    KU{end+1} = ku ;
    KUf{end+1} = ku ; %hrwa.lsfit(S(:,:,ff),ku,fun) ;
    FU{end+1} = 0*ku + ff ;
    
    ester = ESTER(:,ff) ;
    ru = find(ester>=tolEster,1,'last') ;
    KE{end+1} = k{ru}(:,ff) ;
    FE{end+1} = KE{end}*0+ff ;
    
end
KU = cat(1,KU{:}) ; FU =  cat(1,FU{:}) ; KUf =  cat(1,KUf{:}) ;
KE = cat(1,KE{:}) ; FE =  cat(1,FE{:}) ;
plot(real(KU),FU,'.')
plot(real(KE),FE,'.')
set(gca,'zdir','reverse') ;
%%
cla ; %axis equal
%plot3(real(K(:)),f(:)/nF,imag(K(:)),'.') ;
plot3(real(KU(:)),FU(:)/nF,imag(KU(:)),'.') ;
plot3(real(KUf(:)),FU(:)/nF,imag(KUf(:)),'.') ;
plot3(real(KE(:)),FE(:)/nF,imag(KE(:)),'.')%'o','linewidth',.1,'markersize',4) ;
set(gca,'zlim',.05*[-1 0])
end
