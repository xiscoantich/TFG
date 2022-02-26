function freq = FFT(signal)
    %Type of ft: matlab, matrix, dft
    type_ft = 'matrix'; %
    switch type_ft
        case 'matlab'
            freq = fft(signal);
        case 'matrix'
            freq = fft_matrix(signal);
        case 'dft'
            freq = dft(signal);
    end
end

function y = fft_matrix(x)
    x = makepowerof2(x);
    y = FFTCT_matrix_recursive(x);
    if length(x(1,:))~=1 %Check if its a row
        y = reshape(y,1,[]);
    end
end

function y = FFTCT_matrix_recursive(x)
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
        I = eye(m);
        y = [I, I; I, -I]*[zt; zb];
    end
end

function y = makepowerof2(x)
    N = length(x);
    y = x;
    while mod(log(N)/log(2),1)~=0 
        y(N+1) = 0;
        N = N+1;
    end
end

function y = dft(x)
    n = length(x);
    y = NaN(size(x));
    for k = 0:n-1  % For each output element
        s = 0;
        for t = 0:n-1  % For each input element
            s = s+x(t+1)*exp(-2i*pi*t*k/n);
        end
        y(k+1) = s;
    end
end