function [K0,K1,K2,M] = matrices(MESH,MAT)
% BUILD THE SAFE MATRICES ASSOCIATED TO THE 1D Waveguide Porblem
    
% Build the interpolation matrices
    MESH = safe.mesh.interpolation(MESH) ; % generate interpolation matrices

% Kinematic relations: E = [E11;E22;E33;2E12;2E13;2E23] = (B0+k*Bk)*[u1;u2]
    [N,Dx,O] = deal(MESH.N,MESH.Dx,MESH.O) ;
    B0 = [  ...
            Dx{1} O O ; .... E11
            O Dx{2} O ; ... E22
            O O O ; ... E33
            Dx{2} Dx{1} O ; ... 2E12
            O O Dx{1} ; ... 2E13
            O O Dx{2} ; ... 2E23
         ] ;
    Bk = -1i*[  ...
            O O O ; .... E11
            O O O ; ... E22
            O O N ; ... E33
            O O O ; ... 2E12
            N O O ; ... 2E13
            O N O ; ... 2E23
         ] ;
                         
% Stiffness matrices K(k) = K0 + k.K1 + k².K2
    CW = kron(sparse(MAT.C),MESH.Wq) ;
    K0 = B0'*CW*B0 ;
    K1 = Bk'*CW*B0 + B0'*CW*Bk ;
    K2 = Bk'*CW*Bk ;
    
% Mass matrix
    Mi = MAT.rho*N'*MESH.Wq*N ;
    Mi = .5*(Mi+Mi') ; % ensure Symmetry
    M = blkdiag(Mi,Mi,Mi) ; 