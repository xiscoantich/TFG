function [y] = IFFTCT(x)
%  Codigo basado en la demostracion de este link:
%  https://adamsiembida.com/how-to-compute-the-ifft-using-only-the-forward-fft/
n_before_padding = length(x);
x = makepowerof2(x);
N = length(x);
y = (2/N)*conj(FFTCT(conj(x))); %Esta fft no puede tener el corte de frecuencias negativas!
y = y(1:n_before_padding,:);  % get rid of padding before returning
end

function [y] = makepowerof2(x)
N = length(x);
y = x;
while mod(log(N)/log(2),1)~=0 
    y(N+1) = 0;
    N = N+1;
end
end

