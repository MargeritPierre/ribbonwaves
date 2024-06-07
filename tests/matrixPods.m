clc,clear all
q = 100 ;
p = 1 ;
r = 1 ;

S = randn(q,p) + 1i*randn(q,p) ;

% Hankel-Toeplitz block matrix
m = floor(q/3) ; n = q-2*m+1 ;
i1 = (m+1:2*m)'+(0:n-1) ; % hankel
i2 = flip(1:m)'+(0:n-1) ; % toeplitz
X1 = reshape(S(i1,:),m,[]) ; % hankel
X2 = reshape(S(i2,:),m,[]) ; % toeplitz
X = X1 + X2 ;

U = randn(m,r) + 1i*randn(m,r) ;
V = X*X'*U ;

% Implicit products
Fs = fft(S,[],1) ;
Fs = permute(Fs(:,:,:),[1 4 3 2]) ; % [q 1 sz p]

% Hankel-vector product
Wh = X1*X1'*U ;
W = U ;
% W = X1'*W
    W = flip(conj(W),1) ; % [m k o]
    Fw = fft(W,q,1) ; % [q k o]
    Fw = Fw.*Fs ; % [q k o p]
    W = ifft(Fw,[],1) ; % [q k o p]
    W = W(2*m:q,:,:,:) ; % [n k o p]
    W = conj(W) ; % [n k o p]
% W = X1*W ;
    W = flip(W,1) ; % [n k o p]
    Fw = fft(W,q,1) ; % [q k o p]
    Fw = sum(Fw.*Fs,4) ; % [q k o]
    W = ifft(Fw,[],1) ; % [q k o]
    W = W(q-m+1:q,:,:) ; % [m k o]
    
    clf ; plot(real(Wh)) ; plot(real(W),'o')

% Toeplitz-vector product
Wt = X2*X2'*U ;
W = U ;
% W = X2'*W
    W = conj(W) ; % [m k o]
    Fw = fft(W,q,1) ; % [q k o]
    Fw = Fw.*Fs ; % [q k o p]
    W = ifft(Fw,[],1) ; % [q k o p]
    W = W(m:q-m,:,:,:) ; % [n k o p]
    W = conj(W) ; % [n k o p]
% W = X2*W ;
    W = flip(W,1) ; % [n k o p]
    Fw = fft(W,q,1) ; % [q k o p]
    Fw = sum(Fw.*Fs,4) ; % [q k o]
    W = ifft(Fw,[],1) ; % [q k o]
    W = flip(W,1) ;
    W = W(m+1:2*m,:,:) ; % [m k o]
    
    clf ; plot(real(Wt)) ; plot(real(W),':')

    
% Mixed H/T product
Wm = X*X'*U ;
W = U ;
% W1 = X1'*W
    W1 = flip(conj(W),1) ; % [m k o]
    Fw1 = fft(W1,q,1) ; % [q k o]
    Fw1 = Fw1.*Fs ; % [q k o p]
    W1 = ifft(Fw1,[],1) ; % [q k o p]
    W1 = W1(2*m:q,:,:,:) ; % [n k o p]
    W1 = conj(W1) ; % [n k o p]
% W2 = X2'*W
    W2 = conj(W) ; % [m k o]
    Fw2 = fft(W2,q,1) ; % [q k o]
    Fw2 = Fw2.*Fs ; % [q k o p]
    W2 = ifft(Fw2,[],1) ; % [q k o p]
    W2 = W2(m:q-m,:,:,:) ; % [n k o p]
    W2 = conj(W2) ; % [n k o p]
% W = X'*U = W1 + W2
    W = W1 + W2 ;
% W1 = X1*W ;
    W1 = flip(W,1) ; % [n k o p]
    Fw1 = fft(W1,q,1) ; % [q k o p]
    Fw1 = sum(Fw1.*Fs,4) ; % [q k o]
    W1 = ifft(Fw1,[],1) ; % [q k o]
    W1 = W1(q-m+1:q,:,:) ; % [m k o]
% W2 = X2*W ;
    W2 = flip(W,1) ; % [n k o p]
    Fw2 = fft(W2,q,1) ; % [q k o p]
    Fw2 = sum(Fw2.*Fs,4) ; % [q k o]
    W2 = ifft(Fw2,[],1) ; % [q k o]
    W2 = flip(W2,1) ;
    W2 = W2(m+1:2*m,:,:) ; % [m k o]
% W = X*U = W1 + W2
    W = W1 + W2 ;
    
    clf ; plot(real(Wm)) ; plot(real(W),'o')
    %norm(Wm(:)-W(:)) 
    
% Mixed product with common operations
Wm = X*X'*U ;
W = U ;
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
    
    clf ; plot(real(Wm)) ; plot(real(W),'o')
    norm(Wm(:)-W(:)) 
