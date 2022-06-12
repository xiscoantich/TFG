classdef FourierTransformer < handle
    properties (Access = public)
        transtype
        transmethod
        freq
        originalsize
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        function obj = FourierTransformer(cParams)
            obj.init(cParams);
        end
        
        function obj = directTransform(obj,data)
            obj.originalsize = size(data.signal);
            if ndims(data.signal) == 1
                obj.freq = obj.FFT(obj,data.signal);
            elseif ndims(data.signal) == 2
                obj.freq = obj.FFT2(obj,data.signal);
            end
        end
        
        function rec = inverseTransform(obj,data)
            if size(obj.originalsize,2) == 1
                rec = real(obj.IFFT(obj,data));
            elseif size(obj.originalsize,2) == 2
                rec = real(obj.IFFT2(obj,data));
            end
        end
    end
    
     methods (Access = private, Static)
        
        function y = FFT(obj,signal)
            switch obj.transmethod
                case 'matlab'
                    y = fft(signal);
                case 'matrix'
                    y = obj.fft_matrix(signal);
                case 'dft' %Creo que la dft no funciona muy bien :(
                    y = obj.dft(signal);
%                 case 'stft' %No tengo inversa de este metodo
%                     obj.freq = obj.stft(signal);
                otherwise %Caso en el que no este bien o vacio
                    cprintf('err', 'transmethod: matlab, matrix, dft \n');
            end
        end
        
        function y = FFT2 (obj,x)
            switch obj.transmethod
                case 'matlab'
                    y = fft2(x);
                case 'matrix'
                    m=size(x,2); %columnas
                    for i=1:m
                        x1(:,i)=obj.FFT(obj,x(:,i)); %FFT por columnas
                    end
                    n=size(x1,1); %filas
                    for i=1:n
                        y(i,:)=obj.FFT(obj,x1(i,:)); %FFT por filas
                    end
                case 'dct_8by8'
                    x = obj.makedivisibleby8(x);
                    y = obj.dct_8by8(x);
                case 'dct'
                    y = dct2(x);
                otherwise %Caso en el que no este bien o vacio
                    cprintf('err', 'transmethod: matlab, matrix, dct_8by8, dct \n');
            end
        end
        
        function y = IFFT(obj,x)
            %Type of ift: matlab, FFT, Coley-Tukey, idft
            switch obj.transmethod
                case 'matlab'
                    y = ifft(x);
                case 'matrix'
                    y = obj.ifft_FFT(x);
                case 'dft'
                    y = obj.idft(x);
                otherwise %Caso en el que este vacio o error?
                    cprintf('err', 'transmethod: matlab, matrix, dft \n');
            end
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
            switch obj.transmethod
                case 'matlab'
                    y=ifft2(x);
                case 'matrix' %Este metodo no tengo inversa
                    y=ifft2(x);
                case 'dct_8by8'
                    y = obj.idct_8by8(x);
                case 'dct'
                    y = idct2(x);
                otherwise
                    cprintf('err', 'transmethod: matlab, matrix, dct_8by8, dct \n');
            end
        end
        
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.transtype = cParams.transtype;
            obj.transmethod = cParams.transmethod;
            %obj.originalsize = size(data.signal);
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
        
        function y = dct_8by8(~,x)
            %y = blockproc(x,[8 8],@(blkStruct) dct2(blkStruct.data));
            %y = blockproc(x,[8 8],@(blkStruct) obj.dct_8x8(blkStruct.data));
            %imshow(y)
            
            T = dctmtx(8);
            dct = @(block_struct) T * block_struct.data * T';
            y = blockproc(x,[8 8],dct);
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
    end
end