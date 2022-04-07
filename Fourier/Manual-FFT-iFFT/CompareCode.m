close all;
clearvars;
clc;


X=1:4096;      %% important! the size of input needs to be 2^n

tic             % check the calculating time 
Y1 = FFTm(X);
toc

tic
Y2 = fft(X);       % Calculate use built in function
toc

figure
hold on

plot(Y1,'ko','DisplayName','Code by hand')
plot(Y2,'b*','DisplayName','Built in code')
title('FFT')

%% ifft
tic
X1=iFFTm(X);
toc

tic
X2=ifft(X);
toc

figure
hold on

plot(X1,'ko','DisplayName','Code by hand')
plot(X2,'b*','DisplayName','Built in code')
title('iFFT')