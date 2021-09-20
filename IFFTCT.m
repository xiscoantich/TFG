function [x] = IFFTCT(X)
%  Codigo basado en la demostracion de este link:
%  https://adamsiembida.com/how-to-compute-the-ifft-using-only-the-forward-fft/
N=length(X);
x=(1/N)*conj(FFTCT(conj(X)));
end


