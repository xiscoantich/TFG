classdef Compressor < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        frequencyCut
        waveCut
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
            
            wave = obj.data.wave;
            wCut = obj.cutWave(wave);
            
            s.type = 'FREQUENCY';
            s.freq = freqCut;
            s.wave = wCut;
            
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
        

        
    end
    
end