function s = music(W,NFFT,fun) 
% COMPUTES THE MUSIC SPECTRUM from the signal subspace W
% uses the fact that Wo*Wo' = (1-W*W')
% s = sum(abs(fft(Wo,NFFT,1)).^2,2) 
%   = sum(abs(F*Wo).^2,2) 
%   = diag(F*(Wo*Wo')*F') 
%   = diag(F*(1-W*W')*F') 
%   = Fim (1mn - Wmr conj(Wnr)) conj(Fjn) 1ij
%   = Fim conj(Fim) - Fim (Wmr conj(Wnr)) conj(Fin)
%   = m - Fim (Wmr conj(Wnr)) conj(Fin)
%   = m - sum(abs(fft(W,NFFT,1)).^2,2) 
if nargin<2 ; NFFT = 1001 ; end
if nargin<3 ; fun = 'exp' ; end
m = size(W,1) ;

switch fun
    case 'exp'
        s = fft(W,NFFT,1) ;
        n = m ; % squared norm of exponential basis functions
    case 'cos'
        Fc = dct(eye(m),NFFT,1) ;
        n = sum(abs(Fc).^2,2) ; % squared norm of cosine basis functions
        s = dct(W,NFFT,1) ;
end
s = 1-sum(abs(s).^2,2)./n ;