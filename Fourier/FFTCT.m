function [y]= FFTCT(x)
%DIVIDE & CONQUER
%Obtenemos un vector columna
%Cuidado no puedes hacer las traspuesta del vector obtenido porque
%reordena como le da la gana :')
%x = x - mean(x);
n_before_padding = length(x);
x = makepowerof2(x);
N = length(x);
y = FFTCT_recursive(x);
%y = y*(N/n_before_padding);

y = y(1:n_before_padding,:);  % get rid of padding before returning
%El zero padding hace que la señal acabe siendo mas pequeña, lo que tengo
%que hacer es multiplicarla por un factor que todavia no se cual es
end
    
function [y] = makepowerof2(x)
N = length(x);
y = x;
while mod(log(N)/log(2),1)~=0 
    y(N+1) = 0;
    N = N+1;
end
end

function X = FFTCT_recursive(x)
N=length(x);
%Separate the x[N] into even and odd-indexed subsequences
for r=0:(N/2-1)
    %Even
    n=2*r;
    xe(r+1)=x(n+1);
    %Odd
    n=2*r+1;
    xo(r+1)=x(n+1);
end
if N<=2
    X(1,1)=xe+xo;
    X(2,1)=xe-xo;
else
    Xe=FFTCT_recursive(xe);
    Xo=FFTCT_recursive(xo);  
    for k=0:length(xe)-1
        w=exp(-(2*pi*1i/N)*k); 
        X(k+1,1)=Xe(k+1)+w*Xo(k+1);
        X(k+1+(N/2),1)=Xe(k+1)-w*Xo(k+1);
    end
end
end

function y = fft_neg_frec_cut(x)
%Esta funcion se debe aplicar antes de eliminar el zero padding
N = length (x);
y = [x(1:N/2) zeros(N/2)];
end

