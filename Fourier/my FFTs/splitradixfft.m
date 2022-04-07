function X = splitradixfft(x)
% Split-Radix, Decimation-In-Time (DIT), FFT Computation function.
%
% Attention: 1) Length of input signal x should be a power of 2.
%                   2) Only the first iteration of the algorithm is implemented here.

N = length(x);
X = zeros(1,N);

% Generate indices for sequences e[n], f[n], g[n] and h[n].
ind1 = 1:2:N-1;
ind2 = 2:4:N-2;
ind3 = 4:4:N;

g = x(ind1); % N/2 length signal made of samples: x[0], x[2], x[4], x[6],     ..., x[N-2]. 
h = x(ind2);  % N/4 length signal made of samples: x[1], x[5], x[9], x[13]   ..., x[N-3].
i = x(ind3);   % N/4 length signal made of samples: x[3], x[7], x[11], x[15] ..., x[N-1]. 

% Compute their N/2- and N/4-point DFT's via built-in fft:
 G = fft(g);
 H = fft(h);
   I = fft(i);

  % Create the phasor vector:
 k = 0:N/4-1;
 W = exp(-1i*2*pi*k/N);
 
% Apply the formulas.
X(1:N/4)             = G(1:N/4)         +     W.*H +      (W.^3).*I;
X(N/4+1:N/2)     = G(N/4+1:N/2) - 1i*W.*H + 1i*(W.^3).*I;
X(N/2+1:3*N/4) = G(1:N/4)         -      W.*H  -      (W.^3).*I;
X(3*N/4+1:N)    = G(N/4+1:N/2) + 1i*W.*H - 1i*(W.^3).*I;