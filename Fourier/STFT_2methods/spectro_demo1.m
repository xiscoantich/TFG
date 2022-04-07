% Short-Time Fourier Transform - Method I (Spectrogram) Demo.
% This demo shows how the custom function stft1.m and is used to
% to produce the spectrogram of an input signal.

%   Copyright 2020 - 2030, Ilias S. Konsoulas.

%% Workspace Initialization.

clc; clear; close all;

%% Select and load the signal to be analyzed.
% load('chirp','Fs','y'); x = y;
% load('gong', 'Fs','y'); x = y;
% load timit2.asc;  Fs = 8000; x = timit2; 
% load('train','Fs','y'); x = y;
% load('splat','Fs','y'); x = y;
% load('laughter','Fs','y'); x = y;
 [x Fs] = audioread('andean-flute.wav');

%% Play Back the selected audio signal:
soundsc(x,Fs,24);

%% Signal Normalization.
x = x.'/max(abs(x));  

%% STFT Parameters.
L    = length(x);
N    = 512;   % Selected window size.
M    = 450;   % Selected overlap between successive segments in samples.
Nfft = 512;   % Selected number of FFT points.

[t,f,S] = stft1(x,N,M,Nfft,Fs,'hamm');

%% Plot the Spectrogram.
h = figure('Name','STFT - Method I Demo');
colormap('jet');

[T,F] = meshgrid(t,f/1000); % f in KHz.
surface(T,F,10*log10(abs(S.^2) + eps),'EdgeColor','none');

axis tight;
grid on;
title(['Signal Length: ',num2str(L),', Window Length: ', num2str(N),', Overlap: ', num2str(M), ' samples.']);
xlabel('Time (sec)');
ylabel('Frequency (KHz)');
colorbar('Limits',[-80, 40]);
cbar_handle = findobj(h,'tag','Colorbar');
set(get(cbar_handle,'YLabel'),'String','(dB)','Rotation',0);
zlim([-80 40]);