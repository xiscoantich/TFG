function [X]= FFTCT(x)
%DIVIDE & CONQUER
%Obtenemos un vector columna
%Cuidado no puedes hacer las traspuesta del vector obtenido porque
%reordena como le da la gana :')
N1=length(x);
%Check if its power of 2
bool_power2=false;
while bool_power2==false
    [x,bool_power2]=makepowerof2(x,bool_power2);
end
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
    Xe=FFTCT(xe);
    Xo=FFTCT(xo);  
    for k=0:length(xe)-1
        w=exp(-(2*pi*1i/N)*k); 
        X(k+1,1)=Xe(k+1)+w*Xo(k+1);
        X(k+1+(N/2),1)=Xe(k+1)-w*Xo(k+1);
    end
end
end

% elseif lenght(xe)<1
%     Xe=FFTCT(xe);
%     Xo=FFTCT(xo);
% end
% for k=0:length(xe)-1
%     w=e^-(2*pi*1i/N)*k;
%     X(k)=xe
    
function [x,bool_power2]=makepowerof2(x,bool_power2)
N=length(x);
if mod(log(N)/log(2),1)~=0
    x(N+1)=0;
else
    bool_power2=true;
end
end

