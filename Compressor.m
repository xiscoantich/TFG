classdef Compressor < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       keep
       signal
       data
    end
    
    methods (Access = public)
        
        function obj = Compressor(cParams)
            obj.init(cParams)
        end
        
        function cData = computeCompressedSignal(obj)
            fourier.freq = obj.cutFrequency(obj.data.freq);
            fourier.type = obj.data.type.ft;
            fourier.domain = 'FOURIER';
            fourier.dim = obj.data.dim;
            fourier.originalsize = size(obj.data.signal);
            
            wavelet.wave = obj.cutWave(obj.data.wave);
            wavelet.type = obj.data.type.wt;
            wavelet.domain = 'WAVELET';
            wavelet.dim = obj.data.dim;
            wavelet.wave_info = obj.data.wave_info;
            wavelet.motherwave = obj.data.motherwave;
            wavelet.originalsize = size(obj.data.signal);
            
            if obj.data.dim == 1
                s = struct('Fourier',Data(fourier),'Wavelet',Data(wavelet));
            else
                U = obj.data.U;
                S = obj.data.S;
                V = obj.data.V;
                [Ucut,Scut,Vcut] = obj.cutPCA(U,S,V);
                pca.domain = 'PCA';
                pca.dim = obj.data.dim;
                pca.U = Ucut;
                pca.S = Scut;
                pca.V = Vcut;
                s = struct('Fourier',Data(fourier),'Wavelet',Data(wavelet),'PCA',Data(pca));
            end
            
            obj.computeError(s)
            cData = s;
            
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.keep = cParams.keep;
            obj.data = cParams.data;
        end
        
        function  fCut = cutFrequency(obj,freq)
            if obj.keep == 1
                fCut = freq;
            else
                freqSort = sort(abs(freq(:)));
                thresh = freqSort(floor((1-obj.keep)*length(freqSort)));
                ind    = abs(freq)>thresh;
                fCut = freq.*ind;
            end
        end
        
        function  wCut = cutWave(obj,wave)
            if obj.keep == 1
                wCut = wave;
            else
                waveSort = sort(abs(wave(:)));
                thresh = waveSort(floor((1-obj.keep)*length(waveSort)));
                ind    = abs(wave)>thresh;
                wCut = wave.*ind;
            end
        end
        
        function [Ucut,Scut,Vcut] = cutPCA(obj,U,S,V)
            nx = size(U,1);
            ny = size(V,2);
            
            
            %Aqui el cut todavia no esta bien
            %r = round((nx*ny)*(obj.keep/3));
            %r = round(length(S)*(obj.keep/3));
            
            Ucut = U(:,1:ceil(nx*obj.keep));
            Scut = S(1:ceil(nx*obj.keep),1:ceil(ny*obj.keep));
            Vcut = V(:,1:ceil(ny*obj.keep));
        end
        
        function computeError(obj,s)
                original = obj.data.signal ;
                n = length(fieldnames(s));
                [s.Fourier.error,~] = ssim_index(s.Fourier.signal,original);
                [s.Wavelet.error,~] = ssim_index(s.Wavelet.signal,original);
                if n == 3
                     [s.PCA.error,~] = ssim_index(s.PCA.signal, original);
                end

%             original = obj.data.signal ;
%             n = length(fieldnames(s));
%             s.Fourier.error = immse(s.Fourier.signal,original);
%             s.Wavelet.error = immse(s.Wavelet.signal, original);
%             if n == 3
%                  s.PCA.error = immse(s.PCA.signal, original);
%             end
        end
        
        function [mssim, ssim_map] = ssim_index(img1, img2, K, window, L)
            
            if (nargin < 2 || nargin > 5)
                ssim_index = -Inf;
                ssim_map = -Inf;
                return;
            end
            
            if (size(img1) ~= size(img2))
                ssim_index = -Inf;
                ssim_map = -Inf;
                return;
            end
            
            [M N] = size(img1);
            
            if (nargin == 2)
                if ((M < 11) || (N < 11))
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return
                end
                window = fspecial('gaussian', 11, 1.5);	%
                K(1) = 0.01;                             % default settings
                K(2) = 0.03;                             %
                L = 255;                                 %
            end
            
            if (nargin == 3)
                if ((M < 11) || (N < 11))
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return
                end
                window = fspecial('gaussian', 11, 1.5);
                L = 255;
                if (length(K) == 2)
                    if (K(1) < 0 || K(2) < 0)
                        ssim_index = -Inf;
                        ssim_map = -Inf;
                        return;
                    end
                else
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return;
                end
            end
            
            if (nargin == 4)
                [H W] = size(window);
                if ((H*W) < 4 || (H > M) || (W > N))
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return
                end
                L = 255;
                if (length(K) == 2)
                    if (K(1) < 0 || K(2) < 0)
                        ssim_index = -Inf;
                        ssim_map = -Inf;
                        return;
                    end
                else
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return;
                end
            end
            
            if (nargin == 5)
                [H W] = size(window);
                if ((H*W) < 4 || (H > M) || (W > N))
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return
                end
                if (length(K) == 2)
                    if (K(1) < 0 || K(2) < 0)
                        ssim_index = -Inf;
                        ssim_map = -Inf;
                        return;
                    end
                else
                    ssim_index = -Inf;
                    ssim_map = -Inf;
                    return;
                end
            end
            
            C1 = (K(1)*L)^2;
            C2 = (K(2)*L)^2;
            window = window/sum(sum(window));
            img1 = double(img1);
            img2 = double(img2);
            
            ssim_map = filter2(window, img1, 'valid');        % gx
            w1 = filter2(window, img2, 'valid');              % gy
            w2 = ssim_map.*w1;                                % gx*gy
            w2 = 2*w2+C1;                                     % 2*(gx*gy)+C1 = num1
            w1 = (w1-ssim_map).^2+w2;                         % (gy-gx)^2+num1 = den1
            ssim_map = filter2(window, img1.*img2, 'valid');  % g(x*y)
            ssim_map = (2*ssim_map+(C1+C2))-w2;               % 2*g(x*y)+(C1+C2)-num1 = num2
            ssim_map = ssim_map.*w2;                          % num
            img1 = img1.^2;                                   % x^2
            img2 = img2.^2;                                   % y^2
            img1 = img1+img2;                                 % x^2+y^2
            
            if (C1 > 0 && C2 > 0)
                w2 = filter2(window, img1, 'valid');           % g(x^2+y^2)
                w2 = w2-w1+(C1+C2);                            % den2
                w2 = w2.*w1;                                   % den
                ssim_map = ssim_map./w2;                       % num/den = ssim
            else
                w3 = filter2(window, img1, 'valid');           % g(x^2+y^2)
                w3 = w3-w1+(C1+C2);                            % den2
                w4 = ones(size(w1));
                index = (w1.*w3 > 0);
                w4(index) = (ssim_map(index))./(w1(index).*w3(index));
                index = (w1 ~= 0) & (w3 == 0);
                w4(index) = w2(index)./w1(index);
                ssim_map = w4;
            end
            
            mssim = mean2(ssim_map);
            
            return
       
        end
    end
end