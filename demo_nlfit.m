clc
clear all

nX = 100 ;
nY = 10 ;
R = 3 ;
fun = 'exp' ;
SNR = 3e-1 ; 

k = 1*pi*((rand(R,1)-.5) + .0i*(rand(R,1)-.5)) ; 
a = randn(R,nY) + 1i*randn(R,nY) ; ones(R,nY) ;

x = (0:nX-1)' ;
switch fun
    case 'exp'
        s = sum(reshape(a.',[1 nY R]).*exp(1i*reshape(k,[1 1 R]).*x),3) ;
    case 'cos'
        k = abs(real(k))-1i*abs(imag(k)) ;
        s = sum(reshape(a.',[1 nY R]).*cos(reshape(k,[1 1 R]).*x),3) ;
end

noise = randn(size(s)) + 1i*randn(size(s)) ;
s = s + noise./norm(noise,'fro').*norm(s,'fro')./SNR ;

[W,lmbda] = hrwa.subspace(s,R,fun) ;
W = W.*permute(sqrt(lmbda),[2 1 3]) ; % keep the eigenvalue info in the signal subspace

NFFT = 10*nX ;
Fw = hrwa.music(W,NFFT,fun) ;
switch fun
    case 'exp'
        kw = (1-1/NFFT)*pi*linspace(-1,1,NFFT) ;
        Fw = fftshift(Fw,1) ;
    case 'cos'
        kw = (1-1/NFFT)*pi*linspace(0,1,NFFT) ;
end

[ke,Ue,Sre] = hrwa.esprit(W,fun,s) ;

[kf,Uf] = hrwa.lsfit(s,ke,fun) ;
errESPRIT = norm(sort(ke)-sort(k))
errLSFIT = norm(sort(kf)-sort(k))

clf ; plot3(x,real(s),imag(s)) ;

clf 
plot(kw,log10(abs(Fw)))
plot(real(k).'.*[1;1],get(gca,'ylim')') ;
plot(real(ke).'.*[1;1],get(gca,'ylim')',':') ;
plot(real(kf).'.*[1;1],get(gca,'ylim')','--') ;
