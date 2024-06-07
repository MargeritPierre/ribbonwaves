function [wips,wipa,wops,wopa] = cutoffW(GEO,MAT,n_or_wrange)
% Compute the first cutoff frequencies w related to a zero wavenumber k=0
% truncated to the n(integer) first frequencies or in the range wrange(float)
if isinteger(n_or_wrange) ; n = 0:double(n_or_wrange)-1 ; wrange = [-inf inf] ; 
else ; n = 0:1000 ; wrange = n_or_wrange ; end
if numel(wrange)==1 ; wrange = [0 wrange] ; end
n = n(:) ;

[h,b,Q,G,nu,rho,xi] = ribbon.parseInput(GEO,MAT) ;

% IN-PLANE MOTION: 
    % phase velocities
    vl = sqrt(Q./rho) ;
    vt = sqrt(G./rho) ;
    % cutoff frequencies
    [wl,nl] = inrange(@(n)n*pi*vl/b,wrange,n) ;
    [wt,nt] = inrange(@(n)n*pi*vt/b,wrange,n) ;
    wips = sort([wt(mod(nt,2)==0) ; wl(mod(nl,2)==1)],'comparisonmethod','abs') ;
    wipa = sort([wt(mod(nt,2)==1) ; wl(mod(nl,2)==0)],'comparisonmethod','abs') ;

%% OUT-OF-PLANE MOTION
% quantities of interest
    vinf = sqrt(xi).*vt ;
    w0 = sqrt(12)/h*vinf ;
% Transverse shear modes
    [wst,nst] = inrange(@(n)sqrt(w0.^2 + (vt.*n.*pi/b).^2),wrange,n) ; % {sin|cos}(b/2*must)=0
    wops = sort([0 ; wst(mod(nst,2)==0)],'comparisonmethod','abs') ;
    wopa = sort([0 ; wst(mod(nst,2)==1)],'comparisonmethod','abs') ;
    return
    
%% TESTS
    w = 2*pi*fc*linspace(0,5,10000) ;
    k = 0 ;
    % wavenumbers (isotropic part mu_i² = k_i² + k²)
    must2 = (w.^2-w0.^2)./vt.^2 ; % transverse shear
    alpha = (vl./vinf).^2 + 1 ;
    beta2 = (1-(w0./w).^2).*(vl./vinf).^2 ;
    gamma = sqrt(alpha.^2 - 4.*beta2) ;
    mub2 = .5*(w./vl).^2.*(alpha + gamma) ; % bending
    musa2 = .5*(w./vl).^2.*(alpha - gamma) ; % axial shear
    % width-direction wavenumbers & trigo
    kb = sqrt(mub2) ; cb = cos(GEO.b./2.*kb) ; sb = sin(GEO.b./2.*kb) ; 
    ksa = sqrt(musa2) ; csa = cos(GEO.b./2.*ksa) ; ssa = sin(GEO.b./2.*ksa) ;
    kst = sqrt(must2) ; cst = cos(GEO.b./2.*kst) ; sst = sin(GEO.b./2.*kst) ;
    % Components
    Sst1 = -k ;
    Sst2 = kst.^2 - k.^2  ;
    Ast3 = (nu-1).*kst.*k ;
    Sb3 = w0.^2.*(nu.*k.^2 + kb.^2) ;
    Ssa3 = w0.^2.*(nu.*k.^2 + ksa.^2) ;
    Ab1 = kb.*(w.^2-vl.^2.*mub2) ;
    Ab2 = 2*kb.*k.*w0.^2 ;
    Asa1 = ksa.*(w.^2-vl.^2.*musa2) ;
    Asa2 = 2*ksa.*k.*w0.^2 ;
    % precomputation
    Db =  Sb3.*(Asa1.*Sst2-Asa2.*Sst1) ...
         ... kb.*(alpha - 2 - gamma) ...
         ... sqrt(alpha + gamma).*(alpha - 2 - gamma) ...
         ; 
    Dsa =  Ssa3.*(Ab1.*Sst2-Ab2.*Sst1) ...
          ... ksa.*(alpha - 2 + gamma) ...
          ... sqrt(alpha - gamma).*(alpha - 2 + gamma) ...
          ; 
    kbksa = (w.^2./vl.^2).*sqrt(beta2) ;
    kbpksa = w./vl.*sqrt(alpha+2.*sqrt(beta2)) ;
    kbmksa = w./vl.*sqrt(alpha-2.*sqrt(beta2)) ;
    a = vl./vinf ; b = sqrt(1-w0.^2./w.^2) ;
    %beta = 
    DbpDsa = ... Db+Dsa ...
             ... (kb+ksa).*(w.^2-vl.^2.*kb.*ksa) ...
              (kbpksa).*(w.^2-vl.^2.*kbksa) ...
             ... (2-sqrt(beta)).*sqrt(alpha+sqrt(beta)) ...
             ...(kb+ksa).*(kb.^2 + ksa.^2 - alpha.*kb.*ksa) ...
             ... (kb+ksa).*(alpha-2) - sqrt(alpha.^2 - beta).*(kb-ksa) ...
             ... sqrt(alpha+sqrt(beta)).*(alpha-2) - sqrt((alpha.^2 - beta).*(alpha-sqrt(beta))) ...
             ; % Db+Dsa
    DbmDsa = ... Db-Dsa ...
             ... (kb-ksa).*(w.^2+vl.^2.*kb.*ksa) ...
              (kbmksa).*(w.^2+vl.^2.*kbksa) ...
             ... (2+sqrt(beta)).*sqrt(alpha-sqrt(beta)) ...
             ...(kb-ksa).*(kb.^2 + ksa.^2 + alpha.*kb.*ksa) ...
             ... (kb-ksa).*(alpha-2) - sqrt(alpha.^2 - beta).*(kb+ksa) ...
             ... sqrt(alpha-sqrt(beta)).*(alpha-2) - sqrt(alpha.^2 - beta).*sqrt(alpha+sqrt(beta)) ...
             ; % Db-Dsa
    % determinants
%     detBs = (cb.*ssa.*Db - csa.*sb.*Dsa) ; % symmetric modes
%     detBa = (sb.*csa.*Db - ssa.*cb.*Dsa) ; % antisymmetric modes
    detBs = .5*(DbmDsa.*sin(GEO.b./2.*(kbpksa)) - DbpDsa.*sin(GEO.b./2.*(kbmksa))) ; % symmetric modes
    detBa = .5*(DbmDsa.*sin(GEO.b./2.*(kbpksa)) + DbpDsa.*sin(GEO.b./2.*(kbmksa))) ; % antisymmetric modes
    % Analytic cases
    m = ribbon.model(GEO,MAT,0,w) ;
%     clf ;
%     plot(w/wc,real(vl.^2.*kb.*ksa)) ; plot(w/wc,imag(vl.^2.*kb.*ksa)) ;
%     plot(w/wc,real(w.^2)) ; plot(w/wc,imag(w.^2)) ;
%     %%
%     plot(w/wc,real(kbksa),':') ; plot(w/wc,imag(kbksa),':') ;
    clf ;
    plot(w/wc,real(kb+ksa)) ; plot(w/wc,imag(kb+ksa)) ;
    %plot(w/wc,real(kb-ksa)) ; plot(w/wc,imag(kb-ksa)) ;
    plot(w/wc,real(sqrt(w./vl).*sqrt(2i)).*12^.25,':') ; %plot(w/wc,imag(kbpksa),':') ;
    %plot(w/wc,real(kbpksa),':') ; plot(w/wc,imag(kbpksa),':') ;
% %     kk = .5*(w.^2./vl.^2).*(beta)
% %     plot(w/wc,real(kk),'-') ; %plot(w/wc,imag(kk),'-') ;
% %     clf ; 
% %     plot(w/wc,real(kb)) ; plot(w/wc,imag(kb))
% %     plot(w/wc,real(ksa)) ; plot(w/wc,imag(ksa))
% %     plot(w/wc,w./vl/sqrt(2).*real(alpha)) ; plot(w/wc,w./vl/sqrt(2).*imag(alpha))
%     legend
    %%
    clf ;
    plot(w/2/pi/fc,abs(m.Dops./sst./kb./ksa./kst.^2./w0.^2));%./w.^3.*vl*2))
    plot(w/2/pi/fc,abs(m.Dopa./cst./kb./ksa./kst.^2./w0.^2));%./w.^3.*vl*2))
    plot(w/2/pi/fc,abs(detBs),':')
    plot(w/2/pi/fc,abs(detBa),':')
    set(gca,'yscale','log')
    
end

function [f,n] = inrange(fcn,range,n)
    f = fcn(n) ;
    ii = real(f)>=min(real(range)) & real(f)<=max(real(range)) ;
    f = f(ii) ; n = n(ii) ;
end

