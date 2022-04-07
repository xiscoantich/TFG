clc;
close all;
clearvars;

load('chirp','Fs','y'); x3 = y;

load('train','Fs','y'); x4 = y;


fs = Fs;
t = 0:1/fs:2;


x1 = 11*sin(2000*2*pi*t);
x2 = chirp(t,0,1,2000);
figure
plot(t(1:1000),x2(1:1000))
xlim([0 0.1]);
title('Linear Chirp')
xlabel('Time')
ylabel('Amplitude')

figure
tiledlayout(2,2,'TileSpacing','tight')
nexttile
spectrogram(x1,256,250,256,fs,'yaxis')
title('Sine wave')

nexttile
spectrogram(x2,256,250,256,fs,'yaxis')
title('Linear chirp')

nexttile
spectrogram(x3,256,250,256,fs,'yaxis')
title('Bird chirp')

nexttile
spectrogram(x4,256,250,256,fs,'yaxis')
title('Train whistle')




