function d = ger_fft(x) 
% German Style Decimation-In-Time (DIT) radix-2 FFT.
% It comes from pseudocode presented in the following book:
% Robert Plato. Numerische Mathematik kompakt. Vieweg, 4th Ed., 2010, page 53.
% University of Siegen.
  
    N = length(x);
    d = x(bitrevorder(1:N));
    q = log2(N);
    
    for r = 0:q-1        
         M  = 2^r;   th = exp(-pi*1i/M);        
         for k = 1:M            
               for j = 0:2^(q-r-1)-1                                   
                    x = th^(k-1)*d(2*j*M+M+k);
                    d(2*j*M+M+k) = d(2*j*M+k) - x;
                    d(2*j*M+k)      = d(2*j*M+k) + x;                   
               end       
         end
    end