function [t, f, S] = stft2(x,N,M,Nfft,Fs,win_type)
%STFT2 Short-Time Fourier Transform (STFT) - Method II.
%
%   [t, f, S] = stft2(x,N,M,Nfft,Fs,win_type) calculates the Short-Time Fourier Transform 
%   (STFT) or in other words the spectrogram of the input signal x.
%
%   Inputs:
%            x is a row vector that contains the signal to be examined.
%            N is the selected length of the window in samples.
%            M is the selected amount of overlap between successive windows in samples.
%            Nfft is the selected number of DFT calculation points (Nfft>=N). 
%            win_type is a string containing one of the windows shown
%            below. The default value of this variable corresponds to the rectangular window.
%
% Outputs: S is a matrix that contains the complex spectrogram of x.    
%            i.e. S(:,k) = fft(k-th frame) computed on Nfft points.
%          f is a vector of frequency points in Hz where the spectrogram is calculated.
%          t is a vector of time points in sec. Each point corresponds to a
%            signal frame.
%
%   Copyright 2020-2030, Ilias S. Konsoulas.

%% Window selection and contruction.
switch win_type 
    
    case  'cheby'
         win = chebwin(N).';

    case 'blackman'
         win = blackman(N).';

    case 'hamm'
         win = hamming(N).';
        
    case 'hann'
         win = hanning(N).';  
        
    case 'kaiser'
         beta = 5;
         win = kaiser(N,beta).';
         
    case 'gauss'
         win = gausswin(N).';
         
    otherwise  % otherwise use the rectangular window
         win = ones(1,N);
end

%% Input Signal Segmentation Params.
x = x(:).';
L = length(x);

% Number of segments (frames) the signal is divided to.
K = floor((L-M)/(N-M)); 

%% Number of Unique FFT Points.
NUPs = Nfft;   
if isreal(x)
   if mod(Nfft,2)   % if N is odd.
        NUPs = (Nfft+1)/2;
   else             % if N is even.
        NUPs = Nfft/2+1; 
   end
end

%% STFT Calculation
X = zeros(N,K);
S = zeros(Nfft,K);

for k=1:K
    X(:,k) = x((k-1)*(N-M)+1:k*N - (k-1)*M).*win;
    S(:,k) = fft(X(:,k),Nfft);
end

S = S(1:NUPs,:);

%% Frame Time Points
t = (N-1)/Fs + (0:K-1)*(N-M)/Fs; % Frame Time Points.

%% Frequency Points in Hz.
f = (0:NUPs-1)*Fs/Nfft;  


%% NOTES:
% When K is an integer then the following equation holds:
% (N-1)*Ts + (K-1)*(N-M)*Ts = (L-1)*Ts = total signal duration in sec.
% or (N-1) + (K-1)*(N-M) = (L-1).