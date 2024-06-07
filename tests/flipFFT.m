
q = 20 ;
m = 10 ;
r = 1 ;

w = randn(m,1) + 1i*randn(m,1) ;

Fwf = fft(flip(w,1),q,1) ;

Fw = fft(w,q,1) ;

shift = exp(-2i*pi*(0:q-1)'*(q-m+1)./q) ;
Fwfa = Fw.*shift ;
Fwfa = circshift(flip(Fwfa,1),1) ;
norm(diff([Fwfa Fwf],1,2))
% [conj(shift.*Fw) Fwf]
clf ; %plot(real(Fw)) ; 
plot(real(Fwf)) ; plot(imag(Fwf)) ; 
plot(real(Fwfa),':') ; plot(imag(Fwfa),':') ;


%% FFT of a vector flip

% random vector
nx = 24;
x = randn(nx,1) + 1i*randn(nx,1);
fx = fft(x);
% flip and FFT
y = x(end:-1:1);
fy = fft(y);

% multiply with inverse shift operator and multiply
w = exp(-2i*pi*(0:nx-1)'./nx); % inverse cicular shift
fz = circshift(flip(fx.*w),1);
norm(fy-fz)

clf ; 
plot(real(fy)) ; plot(imag(fy)) ; 
plot(real(fz),':') ; plot(imag(fz),':') ;

