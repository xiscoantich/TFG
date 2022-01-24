function [X]= FFT2CT(x)
%FFT2 Cooley-Tukey Algorithm
m=size(x,2); %columnas
for i=1:m
    x1(:,i)=FFTCT(x(:,i)); %FFT por columnas
end
n=size(x1,1); %filas
for i=1:n
    X(i,:)=FFTCT(x1(i,:)); %FFT por filas
end
end

