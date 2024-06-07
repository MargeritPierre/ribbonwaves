function [D,K,W,colors] = detMap(GEO,MAT,k,w,sz,colors)
%% BUILD THE DETERMINANT MAP
if nargin<5 ; sz = [] ; end % optionnal image size
if nargin<6 % colors associated to the different wave groups
    colors = get(0,'DefaultAxesColorOrder') ; % hsv(4)*8/10 + 1/10 ; 
    colors = colors(1:4,:) ;
end

normalize = true ;
filtSz = 5.*[1 1] ;
Dmin = .3 ; Dmax = 1 ; % bounds on the normalized determinant

% Reinterpolate if needed
if ~isempty(sz)
    sz = sz(:)'.*[1 1] ;
    k = interp1([0 cumsum(abs(diff(k)))]/sum(abs(diff(k))),k,linspace(0,1,sz(2))) ;
    w = interp1([0 cumsum(abs(diff(w)))]/sum(abs(diff(w))),w,linspace(0,1,sz(1))) ;
end

% K-W coordinate matrices
if isvector(k) && isvector(w) ; k = k(:).' ; w = w(:) ; end
K = k + w*0 ;
W = k*0 + w ;

% Compute the determinant(s)
m = ribbon.model(GEO,MAT,K,W) ;
Dips = m.Dips_n ; Dipa = m.Dipa_n ; Dops = m.Dops_n ; Dopa = m.Dopa_n ;
                            
% Filter for display
filtFun = @(D)D./medfilt2(abs(D),filtSz,'symmetric') ;
% Concatenate filtered maps
D = cat(4,filtFun(Dips),filtFun(Dipa),filtFun(Dops),filtFun(Dopa)) ;
% Process
D = abs(D) ;
D = max(Dmin,min(D,Dmax)) ; % bound
D = (D-Dmin)./(Dmax-Dmin) ; % rescale

D = 1-sum((1-D).*permute(1-colors,[3 4 2 1]),4) ;

end

