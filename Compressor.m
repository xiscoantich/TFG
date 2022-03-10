classdef Compressor < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        frequencyCut
        waveCut
        Ucut
        Scut
        Vcut
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
            freq = obj.data.freq;             
            freqCut = obj.cutFrequency(freq);
            s.type_ft = obj.data.type_ft;
            
            wave = obj.data.wave;
            wCut = obj.cutWave(wave);
            
            U = obj.data.U;
            S = obj.data.S;
            V = obj.data.V;
            
            [Ucut,Scut,Vcut] = obj.cutPCA(U,S,V);
            
            s.type = 'FREQUENCY';
            s.freq = freqCut;
            s.wave = wCut;
            s.wave_info = obj.data.wave_info;
            s.U = Ucut;
            s.S = Scut;
            s.V = Vcut;
            
            cData = Data(s);
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
            obj.frequencyCut = fCut;
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
            obj.waveCut = wCut;
        end
        
        function [Ucut,Scut,Vcut] = cutPCA(obj,U,S,V)  
            nx = size(U,1);
            ny = size(V,2);
            %r = round((nx*ny)*(obj.keep/3));
            %r = round(length(S)*(obj.keep/3));
            r = 1;
            Ucut = U(:,1:r);
            Scut = S(1:r,1:r);
            Vcut = V(:,1:r);
            obj.Ucut = Ucut;
            obj.Scut = Scut;
            obj.Vcut = Vcut;
        end
        
    end
    
end