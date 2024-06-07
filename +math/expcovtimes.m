function W = expcovtimes(Fs,W)
% product of Cpp = H*H' with W [m k o]
% where H[m n] is a hankel matrix built from the columns of S[q 1 o p]
% Fs is the columnw-wise Fourier transform of S
m = size(W,1) ; q = size(Fs,1) ;
% W = Hp'*W
    W = flip(conj(W),1) ; % [m k o]
    Fw = fft(W,q,1) ; % [q k o]
    Fw = Fw.*Fs ; % [q k o p]
    W = ifft(Fw,[],1) ; % [q k o p]
    W = W(m:q,:,:,:) ; % [n k o p]
    W = conj(W) ; % [n k o p]
% W = Hp*W ;
    W = flip(W,1) ; % [n k o p]
    Fw = fft(W,q,1) ; % [q k o p]
    Fw = sum(Fw.*Fs,4) ; % [q k o]
    W = ifft(Fw,[],1) ; % [q k o]
    W = W(q-m+1:q,:,:) ; % [m k o]
end