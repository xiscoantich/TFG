function y = fft_it2(x)   %#codegen   % Checks whether this function is
                                                             % suitable for automatic .mex code generation.

% This custom function is an iterative implementation of the
% In-Place, Decimation-In-Time (DIT), radix-2 FFT.
% It comes from the paper by Stefan Worner 
% "Fast Fourier Transform. Numerical Analysis Seminar", ETH Zurich.
              
    N= length(x);  
    y = 0; 
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
    
    % This specific script is used to generate the mex function: fft_it2_mex.mexw64.
    % The mex code was generated by the following commands:
    %
    % A) The next 2 commands enable Dynamic Memory Allocation for support
    % of variable-size input data:
    % 
    % >> mexcfg = coder.config('mex');
    % >> mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
    %
    % B) Next Define an "example" input (i.e. a complex row-vector):
    % >> a = randn(1,2^12) + j*randn(1,2^12);
    %
    % C) Finally call the correct codegen version:
    % The expression: typeof(a,[1 Inf]) specifies as an input a row-vector
    % of variable-size to provide maximum functionality:
    %
    % >> codegen -config mexcfg fft_it2.m -args {coder.typeof(a,[1 Inf])}
    
    