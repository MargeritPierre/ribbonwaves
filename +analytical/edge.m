function ve = edge(nu,vt)
% EDGE computes the edge wave velocity
if nargin<2 ; vt = 1 ; end
ve = 0.874+0.054*sqrt(nu).*vt ;
end