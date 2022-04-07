function X = dif_fft(x)
% Radix-2, Decimation in Frequency (DIF) FFT Computation.
%
% Attention: 1) Length of input signal x should be a power of 2.
%                  2) Only the first iteration of this algorithm is implemented here!

N = length(x);

k=1:N/2;
phasor = exp(-1i*2*pi*(k-1)/N);
g = x(k) + x(k + N/2);
h = (x(k) - x(k + N/2)).*phasor;

% Compute their N/2-point DFT:
G = fft(g);
H = fft(h);

X = zeros(1,N); 
X(2*k)     = H(k);  % odd part
X(2*k-1) = G(k);  % even part