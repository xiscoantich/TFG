function [X]= FFTCT(x)
%Separate the x[n] into even and odd-indexed subsequences
N=length(x);
for r=0:(N/2-1)
    %Even
    n=2*r;
    xe(r+1)=x(n+1);
    %Odd
    n=2*r+1;
    xo(r+1)=x(n+1);
end
if length(xe)~=length(xo)
   cprintf('err', 'dimension error \n')
end
if N<=2
    X(1)=xe+xo;
    X(2)=xe-xo;
else
    Xe=FFTCT(xe);
    Xo=FFTCT(xo);  
    for k=0:length(xe)-1
        w=exp(-(2*pi*1i/N)*k);
        X(k+1)=Xe(k+1)+w*Xo(k+1);
        X(k+1+(N/2))=Xe(k+1)-w*Xo(k+1);
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
    

