function X = radix2fft(x)
% Radix-2, Decimation-In-Time (DIT), FFT Computation function.
%
% Attention: 1) Length of input signal x should be a power of 2.
%                  2) Only the first iteration of the algorithm is implemented here.

N = length(x);

% Generate indices for sequences g[n] and h[n].
even_time_ind = 1:2:N-1;
  odd_time_ind = 2:2:N;
  
g = x(even_time_ind); % N/2 point signal made of the even-time-indexed samples of x[n]: x[0], x[2], x[4], ..., x[N-2]. 
h = x(odd_time_ind);   % N/2 point signal made of the odd-time-indexed  samples of x[n]: x[1], x[3], x[5], ..., x[N-1].

% Compute their N/2-point DFT's via built-in fft:
G1 = fft(g);
H1 = fft(h);

% Create the periodic extensions of these 2 DFT's:
H = [H1 H1];
G = [G1 G1];

% Create the phasor vector:
 k = 0:N-1;
 W = exp(-1i*2*pi*k/N);

% Apply the formula.
X = G + W.*H;