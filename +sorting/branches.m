function edg = branches(U,lmbda,tolMAC,tolLmbda)
% Build the branches of a parametric eigenvalue problem
% U[nDOF nModes nP]
% lmbda[nModes nP]
if nargin<3 ; tolMAC = .5 ; end
if nargin<4 ; tolLmbda = 3./size(lmbda(:,:),2) ; end
unique = true ;

% Compute the MACs between subsequent eigenvalue sets
MAC = sorting.MAC(U(:,:,2:end),U(:,:,1:end-1)) ;
% Build first Connectivity
[edg,id] = sorting.adjacency(1-abs(MAC),1-tolMAC,unique) ;
% Eigenvalue distances
dlmbda = sorting.distance(permute(lmbda(:,2:end),[3 1 2]),permute(lmbda(:,1:end-1),[3 1 2]),'range') ;
% Cut too "long" edges 
dL = dlmbda(id) ;
edg(dL>tolLmbda,:) = [] ;
    
end



%%
% Ur = U./sqrt(sum(abs(U).^2,1)) ;
% nQ = floor(min(2*size(lmbda,2),size(U,1)/4)) ;
% [Q,q2] = eigs(Ur(:,:)*Ur(:,:)',nQ,'lm') ;
% q = sqrt(abs(diag(q2))) ;
% [~,nQ] = max(-diff(q)./q(1:end-1)) ; % minimum description length
% Q = Q(:,1:nQ) ;
% Uq = reshape(Q\Ur(:,:),[nQ size(U,2:3)]) ;