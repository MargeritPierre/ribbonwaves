function mac = MAC(U,V) 
% Compute the Modal Assurance Criterion (==normalized scalar product)
% between two bases U[nP nModesU szU(3:end)] and V[nP nModesV szV(3:end)]
% MAC[nModesU nModesV max(szU(3:end),szV(3:end))] ;

% Normalize
U = U./sqrt(sum(abs(U).^2,1)) ;
if nargin<2 || isempty(V) 
    V = U ; % compute the auto-MAC by default
else
    V = V./sqrt(sum(abs(V).^2,1)) ;
end

% Size information
nDims = max(ndims(U),ndims(V)) ;
szU = size(U,1:nDims) ;
szV = size(V,1:nDims) ;
if szV(1)~=szU(1) ; error('the first dimension of U and V should match') ; end
szOK = szU==1 | szV==1 | szV==szU ;
if any(~szOK(3:end)) ; error('the sizes of V and U for dimensions>2 should match or be singleton') ; end 

% Compute MAC
U = permute(U,[2,nDims+1,3:nDims,1]) ;
V = permute(V,[nDims+1,2,3:nDims,1]) ;
szMAC = max(size(U,1:nDims+1),size(V,1:nDims+1)) ;
if prod(szMAC)>4e9/16 % if the cross product is more than 4Gb..
    mac = zeros(szMAC(1:nDims)) ;
    for vv = 1:size(V,2)
        mac(:,vv,:) = sum(conj(U(:,:,:,:,:)).*V(:,vv,:,:,:),nDims+1) ;
    end
else
    mac = sum(conj(U).*V,nDims+1) ;
end

end