function y = fft_it(x)
% This custom function is an iterative implementation of the
% In-Place, Decimation-In-Time (DIT) radix-2 FFT.
% It comes from the paper by Stefan Worner 
% "Fast Fourier Transform. Numerical Analysis Seminar", ETH Zurich.
   
    N= length(x); 
    x = x(bitrevorder(1:N));
    q = log2(N);

    for j = 1:q
        m = 2^(j-1);
        d = exp(-pi*1i/m).^(0:m-1);
        for k = 1:2^(q-j)
              s = 2*(k-1)*m+1;   % start-index
              e = 2*k*m;             % end-index
              r = s + (e-s+1)/2;  % middle-index
              y_top       = x(s:r-1);
              y_bottom = x(r:e);
              z = d .* y_bottom;
              y = [y_top + z, y_top - z];
              x(s:e) = y;
        end
    end