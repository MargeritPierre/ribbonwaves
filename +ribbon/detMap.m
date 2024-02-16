function [D,K,W,colors] = detMap(GEO,MAT,k,w,sz,colors)
%% BUILD THE DETERMINANT MAP
if nargin<5 ; sz = [] ; end % optionnal image size
if nargin<6 ; colors = hsv(4)*8/10 + 1/10 ; end % colors associated to the different wave groups

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
K = k(:).' + w(:)*0 ;
W = k(:).'*0 + w(:) ;

% Compute the determinant(s)
[~,detAs,detAa,detBs,detBa] = ribbon.determinant( ...
                                                    GEO.h ...
                                                    ,GEO.b ...
                                                    ,MAT.Q ...
                                                    ,MAT.G ...
                                                    ,MAT.rho ...
                                                    ,W ...
                                                    ,K ...
                                                    ,normalize ) ;
                            
% Filter for display
filtFun = @(D)D./medfilt2(abs(D),filtSz,'symmetric') ;
% Concatenate filtered maps
D = cat(4,filtFun(detAs),filtFun(detAa),filtFun(detBs),filtFun(detBa)) ;
% Process
D = abs(D) ;
D = max(Dmin,min(D,Dmax)) ; % bound
D = (D-Dmin)./(Dmax-Dmin) ; % rescale

D = 1-sum((1-D).*permute(1-colors,[3 4 2 1]),4) ;

end

