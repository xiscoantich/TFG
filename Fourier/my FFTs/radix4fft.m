function X = radix4fft(x)
% Radix-4, Decimation-In-Time (DIT), FFT Computation function.
%
% Attention: 1) Length of input signal x should be a power of 2.
%                   2) Only the first iteration of the algorithm is implemented here.

N = length(x);

% Generate indices for sequences e[n], f[n], g[n] and h[n].
ind1 = 1:4:N-3;
ind2 = 2:4:N-2;
ind3 = 3:4:N-1;
ind4 = 4:4:N;

e = x(ind1); % N/4 point signal made of samples: x[0], x[4], x[8],   ..., x[N-4]. 
f  = x(ind2);  % N/4 point signal made of samples: x[1], x[5], x[9],   ..., x[N-3].
g = x(ind3); % N/4 point signal made of samples: x[2], x[6], x[10], ..., x[N-2]. 
h = x(ind4); % N/4 point signal made of samples: x[3], x[7], x[11], ..., x[N-1].

% Compute their N/4-point DFT's via built-in fft:
E1 = fft(e);
F1 = fft(f);
G1 = fft(g);
H1 = fft(h);

% Create the periodic extensions of total size N of these 4 DFT's:
E = [E1 E1 E1 E1];
F = [F1 F1  F1 F1];
G = [G1 G1 G1 G1];
H = [H1 H1 H1 H1];

% Create the phasor vector:
 k = 0:N-1;
 W = exp(-1i*2*pi*k/N);

% Apply the formula.
X = E + W.*F + (W.^2).*G + (W.^3).*H;