function MAT = coefficients(MAT) 
% Fills the lamé, E, nu, G, Q coefficients of a material

% Be sure that E, G, Q and nu are defined
if isfield(MAT,'E') && isfield(MAT,'nu')
    MAT.G = MAT.E./2./(1+MAT.nu) ; % shear modulus
    MAT.Q = MAT.E./(1-MAT.nu^2) ; % Plane stress modulus
elseif isfield(MAT,'E') && isfield(MAT,'G')
    MAT.nu = MAT.E./(2*MAT.G)-1 ; % Poisson ratio
    MAT.Q = MAT.E./(1-MAT.nu^2) ; % Plane stress modulus
elseif isfield(MAT,'Q') && isfield(MAT,'G')
    MAT.nu = 1-2*(G./Q) ; % Poisson ratio
    MAT.E = MAT.Q.*(1-MAT.nu.^2) ; % Young Modulus
end
    
    
% Material 3D stiffness (Voigt notation)
    v = MAT.nu ;
    MAT.C = MAT.E./(1+v)./(1-2*v).*[...
                              1-v v v 0 0 0 ; ...
                              v 1-v v 0 0 0 ; ...
                              v v 1-v 0 0 0 ; ...
                              0 0 0 .5-v 0 0 ; ...
                              0 0 0 0 .5-v 0 ; ...
                              0 0 0 0 0 .5-v ; ...
                             ] ;