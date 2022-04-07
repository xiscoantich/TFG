%Fourier zero padding study
close all;
clearvars;

 folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));   

Fs = 4 ;                      %Frecuencia de muestreo de la senyal
time = 20;                      %Duracion de la señal
[x,t] = signal4(Fs,time);
x=x';
x_zp = makepowerof2(x);

f_matlab = reshape(fft(x),[],1);
f_matlab_zp = reshape(fft(x_zp),[],1);
f_matlab_zp = f_matlab_zp(1:length(f_matlab));
f_FFTCT_matrix = FFTCT_matrix(x);
f_FFTCT = FFTCT(x);
%f_FFTm = reshape(FFTm(x),[],1);
figure

subplot(2,3,1);
error1 = abs(f_matlab-f_FFTCT_matrix);
plot(error1); %Aqui podemos ver como aparece un erorr en fourier
title('FFTCTmatrix vs Matlab FFT')
ylabel('Error')

subplot(2,3,2);
error2 = abs(f_matlab-f_FFTCT);
plot(error2); %Aqui podemos ver como aparece un erorr en fourier
title('FFTCT vs Matlab FFT')
ylabel('Error')

% error3 = abs(f_matlab-f_FFTm);
% figure
% plot(error3); %Aqui podemos ver como aparece un erorr en fourier
% title('Error en los coeficientes de Fourier al hacer FFTm')
% ylabel('Error')
subplot(2,3,3);
error4 = abs(f_FFTCT_matrix-f_FFTCT);
plot(error4); %Aqui podemos ver como aparece un erorr en fourier
title('FFTCTmatrix y FFTCT ')
ylabel('Error')


subplot(2,3,4);
error5 = abs(f_FFTCT_matrix-f_matlab_zp);
plot(error5); %Aqui podemos ver como aparece un erorr en fourier
title('FFTCTmatrix y Matlab zp ')
ylabel('Error')



subplot(2,3,5);
plot(real(f_matlab),'x-','DisplayName','Matlab');hold on;
plot(real(f_matlab_zp),'*-','DisplayName','Matlab zp');hold on;
plot(real(f_FFTCT_matrix),'o-','DisplayName','FFTCT_matrix');hold on;
plot(real(f_FFTCT),'+-','DisplayName','FFTCT');hold on;
legend

figure
rec_matlab = ifft(f_matlab);
rec_matlab_zp = IFFTCT(f_matlab_zp);
rec_FFTCT_matrix = IFFTCT(f_FFTCT_matrix);
rec_FFTCT = IFFTCT(f_FFTCT);

plot(real(rec_matlab),'x-','DisplayName','Matlab');hold on;
plot(real(rec_matlab_zp),'*-','DisplayName','Matlab zp');hold on;
plot(real(rec_FFTCT_matrix),'o-','DisplayName','FFTCT_matrix');hold on;
plot(real(rec_FFTCT),'+-','DisplayName','FFTCT');hold on;
legend


function [y] = makepowerof2(x)
    N = length(x);
    y = x;
    while mod(log(N)/log(2),1)~=0
        y(N+1) = 0;
        N = N+1;
    end
end

function [x,t] = signal4(Fs, time) %(Frequencia me muestreo, tiempo de la senyal)
    t = 0:1/Fs:time-(1/Fs);
    t1 = t(1:end/2);            % Primera mitad del vector
    t2 = t(end/2+1:end);        % Segunda mitad del vector
    x=[zeros(1,length(t1)) ones(1,length(t2))];
end