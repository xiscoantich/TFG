function y = IFFT(x,n_before_padding)
%Type of ift: matlab, FFT, Coley-Tukey, idft
    type_ft = 'matrix'; %
    switch type_ft
        case 'matlab'
            y = ifft(x);
        case 'FFT'
            y = ifft_FFT(x);
        case 'Coley-Tukey'
            y = ifft_ct(x);
        case 'idft'
            y = idft(x);
    end
    
    %Cut the zero-padding
    if isempty(n_before_padding)
        return
    else
        if length(x(1,:))~=1
            y = y(:,1:n_before_padding);
        else
            y = y(1:n_before_padding,:);
        end
    end
end

function y = ifft_FFT(x)
%  Codigo basado en la demostracion de este link:
%  https://adamsiembida.com/how-to-compute-the-ifft-using-only-the-forward-fft/
%Este codigo utiliza la funcion FFT, por lo que utiliza el caso que este
%seleccionado

    x = makepowerof2(x);
    N = length(x);
    y = real((1/N)*conj(FFT(conj(x)))); %Esta fft no puede tener el corte de frecuencias negativas!
    % Get rid of padding before returning
end

function y = ifft_ct(x)
    x = makeporof2(x);
    N=max(size(x)); 
    p=log2(N);  
    if p == 1
        y = idft(x);
    return 
    else 
        Even  = iFFTe(x(1:2:end-1));    
        Odd   = iFFTe(x(2:2:end)); 
        Wm   = exp(1i*2*pi/N);
        y    = nan(1,N);
        for i=0:N-1   
            if i<=N/2-1
               y(i+1)= Even(i+1)+Wm^(i)*Odd(i+1);
            else
               y(i+1)= Even(i-N/2+1)+Wm^(i)*Odd(i-N/2+1);
            end        
        end
    end
end

function y = idft(X)    %iDFT  function 
    N = max(size(X));
    Wm = exp(1i*2*pi/N);
    y = nan(1,N);
    int = nan(1,N);
    for k = 0:N-1           
        for x = 0:N-1    
            int(x+1) = Wm^(x*k)*X(x+1);
        end    
        y(k+1) = sum(int);      
    end
end

% function ispowerof2(x)
%     N = length(x);
%     if mod(log(N)/log(2),1)==0
%         return true
%     else
%         return false
%     end
% end

function [y] = makepowerof2(x)
    N = length(x);
    y = x;
    while mod(log(N)/log(2),1)~=0
        y(N+1) = 0;
        N = N+1;
    end
end

