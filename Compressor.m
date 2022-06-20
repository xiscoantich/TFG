classdef Compressor < handle
    
    properties (Access = public)   
        rec
        err
    end
    
    properties (Access = private)
        data
        keep  
        cut
        method
    end
    
    methods (Access = public)
        
        function obj = Compressor(data,cParams)
            obj.init(data,cParams)
        end
        
        function computeErr(obj,odata)
            original = odata.signal ;
            obj.err.mse = immse(original,obj.rec);
            obj.err.mssim_matlab = multissim(obj.rec,original);
            [obj.err.mssim,~] = ssim_index(obj.rec,original);
            obj.err.ssim = ssim(obj.rec,original);
            
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
                
            end
            
        end
    end
    
    methods (Access = private)
        
        function init(obj,data,cParams)
            obj.keep = cParams.keep;
            obj.method = cParams.method;
            obj.data = data;
            obj.computeCompression()
            obj.computeRec()
        end
        
        function computeCompression(obj)
            switch obj.method
                case 'threshold'
                    obj.cut = obj.cutThreshold();
                case 'mask'
                    obj.cut = obj.cutMask();
            end
        end
        
        function cut = cutThreshold(obj)
            
            switch obj.data.transtype
                case 'Fourier'
                    switch obj.data.transmethod
                        case 'dct_8by8'
                            myfun = @(block_struct) Threshold(block_struct.data,obj.keep);
                            block_size = [8 8];
                            cut = blockproc(obj.data.freq,block_size,myfun);
                        otherwise
                            switch obj.data.transmethod
                                case {'matlab', 'matrix'}
                                    keep = obj.keep/2;
                                otherwise
                                    keep = obj.keep;
                            end
                            cut = Threshold(obj.data.freq,keep);
                    end
                case 'Wavelet'
                    cut = Threshold(obj.data.wave,obj.keep);
                case 'PCA'
                    [cut.U,cut.S,cut.V] = cutPCA(obj.data.U,obj.data.S,obj.data.V,obj.keep);
            end
            
            function  cut = Threshold(freq, keep)
                if keep == 1
                    cut = freq;
                else
                    freqSort = sort(abs(freq(:)));
                    thresh = freqSort(floor((1-keep)*length(freqSort)));
                    ind    = abs(freq)>thresh;
                    cut = freq.*ind;
                end
            end
            
            function [Ucut,Scut,Vcut] = cutPCA(U,S,V,keep)
                n = size(U,1);
                m = size(V,1);
                k = keep*(n*m)/(n+m+1);
                Ucut = U(:,1:ceil(k));
                Scut = S(1:ceil(k),1:ceil(k));
                Vcut = V(:,1:ceil(k));
            end
        end
        
        function cut = cutMask(obj)
            switch obj.data.transtype
                case 'Fourier'
                    switch obj.data.transmethod
                        case 'dct_8by8'
                                if obj.keep >= 0.50
                                    n = 8;
                                elseif obj.keep >=0.4375
                                    n = 7;
                                elseif obj.keep >=0.3281
                                    n = 6;
                                elseif obj.keep >= 0.2344
                                    n = 5;
                                elseif obj.keep >=0.1562
                                    n = 4;
                                elseif obj.keep >=0.0938
                                    n = 3;
                                elseif obj.keep >= 0.0469
                                    n = 2;
                                elseif obj.keep >= 1/64
                                    n = 1;        
                                end
                            B = flip(tril(ones(n,n)));
                            mask = zeros (8,8);
                            mask(1:n, 1:n) = B; 
                                
                            cut = blockproc(obj.data.freq,[8 8],@(block_struct) mask .* block_struct.data);
                          
                        otherwise
%                             %No se porque no funciona
%                             %Create circle
%                             imageSizeX= obj.data.originalsize(1);
%                             imageSizeY= obj.data.originalsize(2);
%                             [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
%                             % Next create the circle in the image.
%                             centerX = imageSizeX / 2; % Wherever you want.
%                             centerY = imageSizeY / 2;
%                             npixels = obj.keep*imageSizeX*imageSizeY/2;
%                             radius = round(sqrt(npixels/pi));
%                             circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
%                             cutshift = (fftshift(obj.data.freq))*double(circlePixels);
%                             cut = ifftshift(cutshift);
                            cprintf('err', 'typesignal err:  mask only avaliable for dct_8by8  \n');
                            
                    end
                case 'Wavelet'
                    cprintf('err', 'typesignal err:  mask only avaliable for dct_8by8  \n');
                case 'PCA'
                    cprintf('err', 'typesignal err:  mask only avaliable for dct_8by8  \n');
            end
            
            function circlePixels = createCircle()
                %Create a logical image of a circle with specified
                % diameter, center, and image size.
                % First create the image.
                
                
            end
        end
        
        function computeRec(obj)
            switch obj.data.transtype
                case 'Fourier'
                    ft = FourierTransformer(obj.data);
                    ft.originalsize = obj.data.originalsize;
                    obj.rec = ft.inverseTransform(obj.cut);
                    obj.cutOriginalSize();
                case 'Wavelet'
                    wt = WaveletTransformer(obj.data);
                    wt.originalsize = obj.data.originalsize;
                    obj.rec = wt.inverseTransform(obj.cut);
                    obj.cutOriginalSize();
                    
                case 'PCA'
                    pca = PCATransformer(obj.data);
                    obj.rec = pca.inverseTransform(obj.cut);
            end
        end
        
        function cutOriginalSize(obj)
            if obj.data.originalsize(2) == 1
                n = size(obj.rec);
                n0 = obj.data.originalsize;
                if mod(n0,2)~=0
                    obj.rec = obj.rec(1:end-1);
                    n0 = n0-1;
                end
                if n0 ~= n(1)
                    ncut = floor((n-n0)/2);
                    obj.rec = obj.rec(1+ncut:end-ncut);
                end
            elseif obj.data.originalsize(2) ~= 1
                n = size(obj.rec);
                n0 = obj.data.originalsize;
                if mod(n0(1),2)~=0
                    obj.rec = obj.rec(end-1,:);
                    n0(1) = n0(1)-1;
                end
                if mod(n0(2),2)~=0
                    obj.rec = obj.rec(:,end-1);
                    n0(2) = n0(2)-1;
                end
                if n0(1) ~= n(1)
                    ncut(1) = (n(1)-n0(1))/2;
                    obj.rec = obj.rec(1+ncut(1):end-ncut(1),:);
                end
                if n0(2) ~= n(2)
                    ncut(2) = (n(2)-n0(2))/2;
                    obj.rec = obj.rec(:,1+ncut(2):end-ncut(2));
                end
            end
        end
        
        
       
            
        end
        
            
end