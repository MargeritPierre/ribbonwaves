function W = coscovtimes(Fs,W)
% product of Cpp = X*X' with W [m k o]
% where X[m n]=X1+X2 is a hankel+toeplitz matrix built from the columns of S[q 1 o p]
%     i1 = (m+1:2*m)'+(0:n-1) ; % hankel
%     i2 = flip(1:m)'+(0:n-1) ; % toeplitz
%     X1 = reshape(S(i1,:),m,[]) ; % hankel
%     X2 = reshape(S(i2,:),m,[]) ; % toeplitz
% Fs is the columnw-wise Fourier transform of S
m = size(W,1) ; q = size(Fs,1) ;
% Common
    W = conj(W) ; % [m k o p]
    Fw = fft(W,q,1) ; % [q k o p]
% W1 = X1'*W
    shift = exp(-2i*pi*(0:q-1)'*(q-m+1)./q) ; % [q 1 1]
    Fw1 = Fw.*shift ; % [q k o p]
    Fw1 = circshift(flip(Fw1,1),1) ; % [q k o p]
    shift = exp(2i*pi*(0:q-1)'*m./q) ; % [q 1 1]
    Fw1 = shift.*Fw1 ; % [q k o p]
% Common
    Fw = (Fw+Fw1).*Fs ; % [q k o p]
    W = ifft(Fw,[],1) ; % [q k o p]
    W = W(m:q-m,:,:,:) ; % [n k o p]
    W = conj(W) ; % [n k o p]
    W = flip(W,1) ; % [n k o p]
    Fw = fft(W,q,1) ; % [q k o p]
    Fw = sum(Fw.*Fs,4); % [q k o]
    W = ifft(Fw,[],1) ; % [q k o]
% W1 = X1*W ;
    W1 = W(q-m+1:q,:,:) ; % [m k o]
% W2 = X2*W ;
    W2 = flip(W,1) ;
    W2 = W2(m+1:2*m,:,:) ; % [m k o]
% W = X*U = W1 + W2
    W = W1 + W2 ;
end