function [W,lmbda] = subspace(S,r)
% Vectorized signal subspace estimation via subspace iteration
% S [l p sz], with:
%   - l the number of "data" samples
%   - p the number of data "realizations"
%   - sz other dimensions for vectorization
% Convergence parameters (hard-coded for now..)

% Size info
sz = size(S,1:max(3,ndims(S))) ;
[q,p,sz] = deal(sz(1),sz(2),sz(3:end)) ;
n = floor(q/2)+1 ; m = q-n+1 ; % sizes of the hankel matrix H[M N*p] 
S = permute(S(:,:,:),[1 4 3 2]) ; % [q 1 sz p]

% Convergence
relTolL = 1e-6 ; sqrt(eps) ;
tolW = inf*m*sqrt(eps) ; 
maxIt = 1000 ; 

% Fourier transform of the signal once for all
    Fs = fft(S,[],1) ; % [q 1 sz p]
    
% Use the Lanczos algo
opA = @(U)covtimes(Fs,U) ;
n_k_p = [m r sz] ;
[W,lmbda] = math.eigh(opA,n_k_p,relTolL) ;
return ;


% % Options for the eigenvalue solver
%     opts = [] ;
%     opts.issym = true ;
%     opts.tol = sqrt(eps) ;
% 
% W = NaN(m,r,sz) ;
% lmbda = NaN(1,r,sz) ;
% for xx = 1:prod(sz)
%     X = Fs(:,:,xx,:) ;
%     [W(:,:,xx),lm] = eigs(@(U)covtimes(X,U),m,r,'lm',opts) ;
%     lmbda(:,:,xx) = diag(lm) ;
% end
% return ;

% Subspace order
if nargin<2 ; r = 1 ; end
k = r ; min(2*r,r+8) ; % signal order extension for convergence ?

%if nargin<2 || isempty(W)
    rng("default") % generate repeatable random numbers
    W = randn([m k sz])+1i*randn([m k sz]) ; % random initial subspace
    rng("shuffle") % re-randomize
%end

% Reshape for algo simplicity
    W = W(:,:,:) ; % [m k sz]
    W = qr_mgs(W) ;
%     W = W./sqrt(sum(abs(W).^2,1)) ; % [m k sz]

% Intialize
    lmbda = ones(1,k,size(W,3)) ; % [1 k sz]
    it = 0 ; overtol = true(size(lmbda)) ; % [1 k sz]

% Subspace estimation via subspace iterations W = qr(Cpp*W)
%     clf ; pl = plot(NaN,NaN,'.k') ;
    while any(overtol(:,1:r,:),'all') && it<maxIt
    % Current data to iterate on
        anyover = true(prod(sz),1) ; any(overtol,[1 2]) ; % o = sum(anyover)
        Wi = W(:,:,anyover) ; % subspace vectors to iterate [m k o]
        Fsi = Fs(:,:,anyover,:)  ; % corresponding data [q 1 o p]
    % product of Cpp = Hp*Hp' with W
        Wi = covtimes(Fsi,Wi) ;
    % eigenvalues & normalization
        lmbdai = sqrt(sum(abs(Wi).^2,1)) ; % [1 k o]
        %Wi = Wi./lmbdai ; % [m k o]
    % locking
%         lock = ~overtol(:,:,anyover) ; % [1 k o]
%         lmbdai = lmbdai.*(1-lock) + lock.*lmbda(:,:,anyover) ;
%         Wi = Wi.*(1-lock) + lock.*W(:,:,anyover) ;
    % Sort in descending order
        [lmbdai,is] = sort(lmbdai,2,'descend') ; % [1 k o]
        for oo = 1:size(W,3) 
            Wi(:,:,oo) = Wi(:,is(:,:,oo),oo) ; % [m k o]
%             lock(:,:,oo) = lock(:,is(:,:,oo),oo) ; % [1 k o]
        end
    % QR factorization
        Wi = qr_mgs(Wi) ;
    % Update
        dW = Wi-W(:,:,anyover) ; % [m k o]
        W(:,:,anyover) = Wi ; % [m k sz]
        dL = (lmbdai-lmbda(:,:,anyover))./lmbdai ; % [1 k o]
        lmbda(:,:,anyover) = lmbdai ; % [1 k sz]
    % Convergence tests
        maxdW = max(abs(dW),[],1) ; % [1 r o]
        maxdL = abs(dL) ; % [1 r o]
        overtol(:,:,anyover) = maxdW>tolW | abs(dL)>relTolL ;
        it = it+1 ;
        disp("SUBSPACE:"...
                + " it=" + string(it) ...
                + ", nOverTol=" + string(sum(overtol(:,1:r,:),'all')) ...
                + ", max(dW)=" + string(max(maxdW(:,1:r,:),[],'all')) ...
                + ", max(dL)=" + string(max(maxdL(:,1:r,:),[],'all')) ...
                ) ;
%     pl.XData = [pl.XData(:) ; it + lmbda(:)*0] ; pl.YData = [pl.YData(:) ; lmbda(:)] ;  drawnow ;
    end
    
%         disp("SUBSPACE:"...
%                 + " it=" + string(it) ...
%                 + ", nOverTol=" + string(sum(overtol(:,1:r,:),'all')) ...
%                 + ", max(dW)=" + string(max(maxdW(:,1:r,:),[],'all')) ...
%                 + ", max(dL)=" + string(max(maxdL(:,1:r,:),[],'all')) ...
%                 ) ;

    lmbda = lmbda(:,1:r,:) ;
    W = qr_mgs(W(:,1:r,:)) ;

end

function W = covtimes(Fs,W)
% product of Cpp = H*H' with W [m k o]
% where H[m n] is a hankel matrix built from the columns of S[q 1 o p]
% Fs is the columnw-wise Fourier transform of S
m = size(W,1) ; q = size(Fs,1) ;
% W = Hp'*W
    W = flip(conj(W),1) ; % [m k o]
    Fw = fft(W,q,1) ; % [q k o]
    Fw = Fw.*Fs ; % [q k o p]
    W = ifft(Fw,[],1) ; % [q k o p]
    W = W(m:q,:,:,:) ; % [n k o p]
    W = conj(W) ; % [n k o p]
% W = Hp*W ;
    W = flip(W,1) ; % [n k o p]
    Fw = fft(W,q,1) ; % [q k o p]
    Fw = sum(Fw.*Fs,4) ; % [q k o]
    W = ifft(Fw,[],1) ; % [q k o]
    W = W(q-m+1:q,:,:) ; % [m k o]
end

function W = qr_mgs(W)
% QR factorization of W [m k o] via modified Gram-Schmidt
% Normalize the first vector
    W(:,1,:) = W(:,1,:)./sqrt(sum(abs(W(:,1,:)).^2,1)) ;
% Orthogonalize the remaining vectors
    for ii = 2:size(W,2)
    % reference vector against which to orthogonalize
        wi = W(:,ii-1,:) ; % [m 1 o]
    % remaining vectors
        Wi = W(:,ii:end,:) ; % [m k-ii+1 o]
    % remove projection
        Wi = Wi - sum(conj(wi).*Wi,1).*wi ; % [m k-ii+1 o]
    % normalize
        Wi(:,1,:) = Wi(:,1,:)./sqrt(sum(abs(Wi(:,1,:)).^2,1)) ; % [m 1 o]
    % save
        W(:,ii:end,:) = Wi ; % [m k-ii+1 o]
    end
end






