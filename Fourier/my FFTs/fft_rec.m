function y = fft_rec(x)
% This custom function is a recursive implementation of the
% Decimation-In-Time (DIT), radix-2 FFT.
% It comes from the paper by Stefan Worner: 
% "Fast Fourier Transform. Numerical Analysis Seminar", ETH Zurich.

N = length(x);
phasor = exp(-2*pi*1i/N) .^ (0:N/2-1);

if N == 1
    y = x;
else
    y_top       = fft_rec(x(1:2:(N-1)));
    y_bottom = fft_rec(x(2:2:N));
    z = phasor.*y_bottom;
    y = [ y_top + z , y_top - z ];
end