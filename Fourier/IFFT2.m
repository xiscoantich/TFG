function y = IFFT2 (x)
    n=size(x,1); %filas
    for i=1:n
        x1(i,:)=IFFT(x(i,:)); %FFT por filas
    end
        m=size(x1,2); %columnas
    for i=1:m
        y(:,i)=IFFT(x1(:,i)); %FFT por columnas
    end
end