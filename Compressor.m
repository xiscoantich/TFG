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
    end
end