function [y] = FFTCT_matrix(x)
%The FFT Via Matrix Factorizations
%Charles Van Loan - Department of Computer Science - Cornell University
%x = x - mean(x);
n_before_padding = length(x);
x = makepowerof2(x);
N = length(x);
y = FFTCT_matrix_recursive(x);
%y = fft_neg_frec_cut (y);
y = y(1:n_before_padding,:);  % get rid of padding before returning 
end

function [y] = FFTCT_matrix_recursive(x)
n = length(x);
if n == 1
    y = x;
else
    m = n/2;
    w = exp(-2*pi*1i/n);
    sigma = zeros (m, 1);
    for j=0:1:m-1
        sigma(j+1,1) = w^j; 
    end
    zt = FFTCT_matrix_recursive(x(1:2:n));
    zb = sigma.*FFTCT_matrix_recursive(x(2:2:n));
    
    I = eye(m); %Matriz identidad
    y = [I, I; I, -I]*[zt; zb];
end
end

function y = fft_neg_frec_cut(x)
%Esta funcion se debe aplicar antes de eliminar el zero padding
N = length (x);
y = [x(1:N/2); zeros(N/2,1)];
end

function [y] = makepowerof2(x)
N = length(x);
y = x;
while mod(log(N)/log(2),1)~=0 
    y(N+1) = 0;
    N = N+1;
end
end