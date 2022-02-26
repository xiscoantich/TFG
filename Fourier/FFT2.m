function y = FFT2 (x)
    m=size(x,2); %columnas
    for i=1:m
        x1(:,i)=FFT(x(:,i)); %FFT por columnas
    end
    n=size(x1,1); %filas
    for i=1:n
        y(i,:)=FFT(x1(i,:)); %FFT por filas
    end
end