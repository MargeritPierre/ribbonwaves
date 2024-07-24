classdef polyexp < math.poly
% POLYEXP polynomial expansion reprezentation: same as polynomial but with systematic order truncation
    
methods
    function this = polyexp(N)
    % POLYEXP Construct the polynomial expansion
        this = this@math.poly([0;1;zeros(N-2,1)]) ;
    end
end

end

