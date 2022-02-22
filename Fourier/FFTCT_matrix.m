function [y] = FFTCT_matrix(x)
%The FFT Via Matrix Factorizations
%Charles Van Loan - Department of Computer Science - Cornell University
%x = x - mean(x);
n_before_padding = length(x);
x = makepowerof2(x);
%N = length(x);
y = FFTCT_matrix_recursive(x);

%Metodo 1
y = fft_neg_frec_cut (y); %Esta linea no se si se deberia ocultar y se deberia arreglar porque ni la priema ni la ultima frequencia tienen un dupicado!!
y = y(1:n_before_padding,:);  % get rid of padding before returning 

%Medodo 2
%y = fft_cut_padding(y,n_before_padding);
%Este metoodo dos no funciona, no es la manera correcta

%Extract the result in the same shape as entered 
if length(x(1,:))~=1 %Check if its a row
    y = reshape(y,1,[]);
end
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
%y = [x(1:N/2+1); zeros(N/2+1,1)];

%Invertir las frecuencias positivas para obtener las negativas
%Importante, la DC y la de Nyquist no estan repetidas
y_neg = x(2:N/2);
y_neg_fliped = flipud(conj(y_neg));

y = [x(1:N/2+1); y_neg_fliped];
end

function x_cut = fft_cut_padding(x,n_before_padding)
%Intercala la part positiva y negativa del vector x
%Alternate
N = length (x);
p = x(1:N/2); %Positive
n = x(N/2+1:end); %Negative
n = flipud(n);

x_alternate = zeros(N,1);
x_alternate(1:2:end) = p;
x_alternate(2:2:end) = n;

%Cut
x_alternate_cut = x_alternate(1:n_before_padding,:);

%Return to the original shape
p_cut = x_alternate_cut(1:2:end);
n_cut = x_alternate_cut(2:2:end);
n_cut = flipud(n_cut);
x_cut = [p_cut; n_cut];
end

function [y] = makepowerof2(x)
N = length(x);
y = x;
while mod(log(N)/log(2),1)~=0 
    y(N+1) = 0;
    N = N+1;
end
end