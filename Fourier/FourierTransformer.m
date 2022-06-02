classdef FourierTransformer < handle

    properties (Access = public)
        win_type_stft
        type_ft
    end
    
    methods (Access = public)

        function freq = directTransform(obj,cParams)
            data = cParams.data;
            signal = data.signal;
            obj.type_ft = data.type.ft;
            if data.dim == 1
                freq = obj.FFT(signal);
            else
                freq = obj.FFT2(signal);
            end 
        end
        
        function signal = inverseTransform(obj,cParams)
                data = cParams.data;
                freq = data.freq;
                obj.type_ft = data.type.ft;
                if data.dim == 1
                    signal = obj.IFFT(freq);
                else
                    signal = obj.IFFT2(freq); 
                end
        end
    end
    
    methods (Access = private)
        function freq = FFT(obj,signal)
            %Type of ft: matlab, matrix, dft,stft
            switch obj.type_ft
                case 'matlab'
                    freq = fft(signal);
                case 'matrix'
                    freq = obj.fft_matrix(signal);
                case 'dft' %Creo que la dft no funciona muy bien :(
                    freq = obj.dft(signal);
                case 'stft'
                    freq = obj.stft(signal);
                otherwise %Caso en el que no este bien o vacio??
                    cprintf('err', 'Type of ft: matlab, matrix, dft,stft \n Default: matlab \n');
                    freq = fft(signal); 
            end
        end
        
        function y = fft_matrix(obj,x)
            x = obj.makepowerof2(x);
            y = obj.FFTCT_matrix_recursive(x);
            if length(x(1,:))~=1 %Check if its a row
                y = reshape(y,1,[]);
            end
        end
        
        function y = FFTCT_matrix_recursive(obj,x)
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
                zt = obj.FFTCT_matrix_recursive(x(1:2:n));
                zb = sigma.*obj.FFTCT_matrix_recursive(x(2:2:n));
                I = eye(m);
                y = [I, I; I, -I]*[zt; zb];
            end
        end
        
        function y = makepowerof2(obj,x)
            %Old power of 2
%             N = length(x);
%             y = x;
%             
%             while mod(log(N)/log(2),1)~=0
%             y(N+1) = 0;
%             N = N+1;
%             end
            
            %New power of two
            N = length(x);
            y = x;
            if mod(N,2) ~= 0
                y(end+1)=y(end);
                N=N+1;
            end
              p = nextpow2(N);
              extension = (2^p-N)/2;
              y = wextend('1D','sym',y,extension);
        end
        
        function y = makedivisibleby8 (obj,x)
            N = size(x);
            
            if mod(N(1),2) ~= 0
                y(end+1,:)=y(end,:);
                N(1)=N(1)+1;
            end
            if mod(N(2),2) ~= 0
                y(:,end+1)=y(:,end);
                N(2)=N(2)+1;
            end
            
            b1=0;
            b2=0;
            if rem(N(1),8)~= 0
            b1 = N(1) + (8 - rem(N(1),8));
            b1 = b1-N(1);
            end
            if rem(N(2),8) ~= 0
            b2 = N(2) + (8 - rem(N(2),8));
            b2 = b2-N(2);
            end
            x1 = wextend('addrow','zpd',x,b1/2);
            y = wextend('addcol','zpd',x1,b2/2);
            
        end
        
        function y = dft(obj,x)
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
        
        function S = stft (obj,x)
            
%            Inputs:
%            x is a row vector that contains the signal to be examined.
%            N is the selected length of the window in samples.
%            M is the selected amount of overlap between successive windows in samples.
%            Nfft is the selected number of DFT calculation points (Nfft>=N). 
%            win_type is a string containing one of the windows shown
%            below. The default value of this variable corresponds to the rectangular window.
            
            %Default parameters???
            N    = 256;   % Selected window size.
            M    = 220;   % Selected overlap between successive segments in samples.
            Nfft = 512;   % Selected number of FFT points.

            if isempty(obj.win_type_stft)
                obj.win_type_stft = 'rectangular';
            end
            switch obj.win_type_stft
                case  'cheby'
                    win = chebwin(N).';
                    
                case 'blackman'
                    win = blackman(N).';
                    
                case 'hamm'
                    win = hamming(N).';
                    
                case 'hann'
                    win = hanning(N).';
                    
                case 'kaiser'
                    beta = 5;
                    win = kaiser(N,beta).';
                    
                case 'gauss'
                    win = gausswin(N).';
                    
                otherwise  % otherwise use the rectangular window
                    win = ones(1,N);
            end
            
            % Input Signal Segmentation Params.
            x = x(:).';
            L = length(x);
            
            % Number of segments (frames) the signal is divided to.
            K = floor((L-M)/(N-M));
            
            % Number of Unique FFT Points.
            NUPs = Nfft;
            if isreal(x)
                if mod(Nfft,2)   % if N is odd.
                    NUPs = (Nfft+1)/2;
                else             % if N is even.
                    NUPs = Nfft/2+1;
                end
            end
            
            % STFT Calculation
            X = zeros(N,K);
            S = zeros(Nfft,K);
            
            for k=1:K
                X(:,k) = x((k-1)*(N-M)+1:k*N - (k-1)*M).*win;
                S(:,k) = fft(X(:,k),Nfft);
            end
            
            S = S(1:NUPs,:);
            
            obj.plotSpectogram(S,L,N,M);
        end
        
        function y = dct_8by8(obj,x)
            y = blockproc(x,[8 8],@(blkStruct) dct2(blkStruct.data));
            %y = blockproc(x,[8 8],@(blkStruct) obj.dct_8x8(blkStruct.data));
            imshow(y)
            
%             T = dctmtx(8);
%             dct = @(block_struct) T * block_struct.data * T';
%             y = blockproc(x,[8 8],dct);
%             imshow(y)
        end
        
        function y = idct_8by8(obj,x)
            y = blockproc(x,[8 8],@(blkStruct) idct2(blkStruct.data));
            
%             mask = [1   1   1   1   0   0   0   0
%                 1   1   1   0   0   0   0   0
%                 1   1   0   0   0   0   0   0
%                 1   0   0   0   0   0   0   0
%                 0   0   0   0   0   0   0   0
%                 0   0   0   0   0   0   0   0
%                 0   0   0   0   0   0   0   0
%                 0   0   0   0   0   0   0   0];
%             B2 = blockproc(x,[8 8],@(block_struct) mask .* block_struct.data);
%             T = dctmtx(8);
%             invdct = @(block_struct) T' * block_struct.data * T;
%             y = blockproc(B2,[8 8],invdct);
%             imshow(mat2gray(y));
        end
        
        function O = dct_8x8(obj,I) %No esta implementada
            cosines = [1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
                0.9808  0.8315  0.5556  0.1951 -0.1951 -0.5556 -0.8315 -0.9808
                0.9239  0.3827 -0.3827 -0.9239 -0.9239 -0.3827  0.3827  0.9239
                0.8315 -0.1951 -0.9808 -0.5556  0.5556  0.9808  0.1951 -0.8315
                0.7071 -0.7071 -0.7071  0.7071  0.7071 -0.7071 -0.7071  0.7071
                0.5556 -0.9808  0.1951  0.8315 -0.8315 -0.1951  0.9808 -0.5556
                0.3827 -0.9239  0.9239 -0.3827 -0.3827  0.9239 -0.9239  0.3827
                0.1951 -0.5556  0.8315 -0.9808  0.9808 -0.8315  0.5556 -0.1951];
            alpha = [0.1250  0.1768  0.1768  0.1768  0.1768  0.1768  0.1768  0.1768
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500];
            O = double(zeros(8,8));
            for p = 1 : 8
                for q = 1 : 8
                    s = double(0);
                    for m = 1 : 8
                        for n = 1 : 8
                            s = s + (double(I(m,n)) * cosines(p,m) * cosines(q,n));
                        end
                    end
                    O(p,q) = alpha(p,q) * s;
                end
            end
        end
        
        function O = idct_8x8(obj,I) %No esta implementada
            cosines = [1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
                0.9808  0.8315  0.5556  0.1951 -0.1951 -0.5556 -0.8315 -0.9808
                0.9239  0.3827 -0.3827 -0.9239 -0.9239 -0.3827  0.3827  0.9239
                0.8315 -0.1951 -0.9808 -0.5556  0.5556  0.9808  0.1951 -0.8315
                0.7071 -0.7071 -0.7071  0.7071  0.7071 -0.7071 -0.7071  0.7071
                0.5556 -0.9808  0.1951  0.8315 -0.8315 -0.1951  0.9808 -0.5556
                0.3827 -0.9239  0.9239 -0.3827 -0.3827  0.9239 -0.9239  0.3827
                0.1951 -0.5556  0.8315 -0.9808  0.9808 -0.8315  0.5556 -0.1951];
            alpha = [0.1250  0.1768  0.1768  0.1768  0.1768  0.1768  0.1768  0.1768
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500
                0.1768  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500  0.2500];
            O = double(zeros(8,8));
            for m = 1 : 8
                for n = 1 : 8
                    s = double(0);
                    for p = 1 : 8
                        for q = 1 : 8
                            s = s + (alpha(p,q) * double(I(p,q)) * cosines(p,m) * cosines(q,n));
                        end
                    end
                    O(m,n) = s;
                end
            end
        end
        
        function y = FFT2 (obj,x)
            switch obj.type_ft
                case 'matlab'
                    y=fft2(x);
                case 'matrix'
                    m=size(x,2); %columnas
                    for i=1:m
                        x1(:,i)=obj.FFT(x(:,i)); %FFT por columnas
                    end
                    n=size(x1,1); %filas
                    for i=1:n
                        y(i,:)=obj.FFT(x1(i,:)); %FFT por filas
                    end
                case 'dct_8by8'
                    x = obj.makedivisibleby8(x);
                    y = obj.dct_8by8(x);
                case 'dct'
                    y = dct2(x);
            end
        end
            
        function y = IFFT(obj,x)
            %Type of ift: matlab, FFT, Coley-Tukey, idft
            switch obj.type_ft
                case 'matlab'
                    y = ifft(x);
                case 'matrix'
                    y = obj.ifft_FFT(x);
                case 'Coley-Tukey'
                    y = obj.ifft_ct(x);
                case 'dft'
                    y = obj.idft(x);
                otherwise %Caso en el que este vacio o error?
                    y = ifft(x);
            end
            
%             %Cut the zero-padding
%             if isempty(n_before_padding)
%                 return
%             else
%                 if length(x(1,:))~=1
%                     y = y(:,1:n_before_padding);
%                 else
%                     y = y(1:n_before_padding,:);
%                 end
%             end
        end
        
        function y = ifft_FFT(obj,x)
            %  Codigo basado en la demostracion de este link:
            %  https://adamsiembida.com/how-to-compute-the-ifft-using-only-the-forward-fft/
            %Este codigo utiliza la funcion FFT, por lo que utiliza el caso que este
            %seleccionado
            x = obj.makepowerof2(x);
            N = length(x);
            y = real((1/N)*conj(obj.FFT(conj(x))));
        end
        
        function y = ifft_ct(obj,x)
            x = makeporof2(x);
            N=max(size(x));
            p=log2(N);
            if p == 1
                y = idft(x);
                return
            else
                Even  = obj.ifft_ct(x(1:2:end-1));
                Odd   = obj.ifft_ct(x(2:2:end));
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
        
        function y = idft(obj,X)    %iDFT  function
            N = max(size(X));
            Wm = exp(1i*2*pi/N);
            y = nan(1,N);
            int = nan(1,N);
            for k = 0:N-1
                for x = 0:N-1
                    int(x+1) = Wm^(x*k)*X(x+1);
                end
                y(k+1) = real(sum(int));
            end
            y = y/length(y);
        end
        
        function y = IFFT2(obj,x)
              %No funcionaa
              
%             n=size(x,1); %filas
%             for i=1:n
%                 x1(i,:)=obj.IFFT(x(i,:)); %FFT por filas
%             end
%             m=size(x1,2); %columnas
%             for i=1:m
%                 y(:,i)=obj.IFFT(x1(:,i)); %FFT por columnas
%             end

%             n=size(x,2); %columnas
%             for i=1:n
%                 x1(:,i)=obj.IFFT(x(:,i)); %FFT por columnas
%             end
%             m=size(x1,1); %filas
%             for i=1:m
%                 y(i,:)=obj.IFFT(x1(i,:)); %FFT por columnas
%             end
            switch obj.type_ft
                case 'dct'
                    y = idct2(x);
                    
                case 'dct_8by8'
                    y = obj.idct_8by8(x);
                    
                otherwise
                    y = ifft2(x);
                    
                
            end
        end
        
    end
    
    methods (Access = private, Static)
        function plotSpectogram(S,L,N,M)
            h = figure('Name','Spectogram');
            colormap('jet');
            
            %[T,F] = meshgrid(t,f/1000); % f in KHz. No tenemos FS
            surface(10*log10(abs(S.^2) + eps),'EdgeColor','none');
            
            axis tight;
            grid on;
            title(['Signal Length: ',num2str(L),', Window Length: ', num2str(N),', Overlap: ', num2str(M), ' samples.']);
            xlabel('Time (sec)');
            ylabel('Frequency'); %Esto deberia ser KHz...
            colorbar('Limits',[-80, 40]);
            cbar_handle = findobj(h,'tag','Colorbar');
            set(get(cbar_handle,'YLabel'),'String','(dB)','Rotation',0);
            %zlim([-80 40]);
        end
    end
end