function C = vmtimes(A,B,mode)
% Vectorized mtimes between matrices A[m,n,p] and X[n,k,p]
% Size info
    nDims = max(ndims(A),ndims(B))+1 ;
    szA = size(A,1:nDims) ; [m,n,szA] = deal(szA(1),szA(2),szA(3:end)) ;
    szB = size(B,1:nDims) ; [n,k,szB] = deal(szB(1),szB(2),szB(3:end)) ;
    sz = max(szA,szB) ;
    p = prod(sz) ;
    
% Usual non-vectorized product
    if p==1 ; C = A*B ; return ; end
% Vectorized versions
    if nargin<3 ; mode = 'loop' ; end
    switch mode
        case 'loop' % Using a loop
        % Initialize
            C = NaN([m k sz]) ;
        % Size indices
            if p==sz(1)
                pA = 1:p ; pB = pA ;
            else
                ii = ind2sub(sz,(1:p)') ;
                iA = min(ii,szA) ;
                iB = min(ii,szB) ;
                pA = sub2ind(szA,iA) ;
                pB = sub2ind(szB,iB) ;
            end
        % Loop
            for pp = 1:p ; C(:,:,pp) = A(:,:,pA(pp))*B(:,:,pB(pp)) ; end
        case 'broadcast' % Using broadcast
            C = reshape(sum(reshape(A,[m n 1 p]).*reshape(B,[1 n k p]),2),[m k sz]) ;
        otherwise ; error('unknown vectorization mode') ;
    end
end

function test
%% UNIT TEST
clc,clear all
m = 10 ;
n = 10 ;
k = 10 ;
szA = [100 1 10] ;
szB = [1 10 10] ;

A = randn([m,n,szA]) ;
B = randn([n,k,szB]) ; 

profile on
tic ;
C = math.vmtimes(A,B) ;
toc
profile off
%profile viewer


%% CHOOSE THE MODE
clc,clear all

m = logspace(1,3,10) ; n = m ; k = 10 ;  p = 10 ; % square mtx
%m = 10 ; n = 1 ; k = 10 ;  p = round(logspace(0,5,10)) ; % Long list of small mtx-vc products

M = round(m) + 0*(n+k+p) ;
N = round(n) + 0*(m+k+p) ;
K = round(k) + 0*(n+m+p) ;
P = round(p) + 0*(n+k+m) ;

Tl = NaN(size(M)) ;
Tb = NaN(size(M)) ;
Tp = NaN(size(M)) ;
for ii=1:numel(M)
    A = randn([M(ii),N(ii),P(ii)]) ;
    B = randn([N(ii),K(ii),P(ii)]) ; 
    disp("VMTIMES | m="+string(M(ii))+",n="+string(N(ii))+",k="+string(K(ii))+",p="+string(P(ii)))
    Tl(ii) = timeit(@()math.vmtimes(A,B,'loop')) ;
    Tb(ii) = timeit(@()math.vmtimes(A,B,'broadcast')) ;
    Tp(ii) = timeit(@()pagemtimes(A,B)) ;
end

clf ; X = M ;
plot(X(:),Tl(:),'o','displayname','LOOP') ; 
plot(X(:),Tb(:),'+','displayname','BROADCAST') ; 
plot(X(:),Tp(:),'x','displayname','PAGEMTIMES') ; 
plot(X(:),1e-9*M(:).*N(:).*K(:).*P(:),':','displayname','FLOPS') ;
set(gca,'xscale','log','yscale','log') ;
legend

end
