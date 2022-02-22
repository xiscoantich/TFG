classdef Data < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        signal
        dim
        freq
    end
    
    properties (Access = private)
        name
    end
    
    methods (Access = public)

        function obj = Data(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            switch cParams.type 
                case 'TEMPORAL'
                   obj.name = cParams.name;
                   obj.loadAudioSignal();
                   obj.dim = size(obj.signal,2);                  
                   obj.computeFourierRepresentation();            

                case 'FREQUENCY'
                   obj.freq = cParams.freq;
                   obj.dim = size(obj.freq,2);                   
                   obj.computeTimeRepresentation();
                   

                case 'WAVELET'

            end
        end

        function  loadImage(obj)
           A = imread(obj.name);
           obj.signal=rgb2gray(A);
        end
        
        function loadAudioSignal(obj)
            if (strcmp(obj.name,'chirp'))
                load('chirp','Fs','y'); 
            elseif (strcmp(obj.name,'gong'))
                load('gong','Fs','y'); 
            elseif (strcmp(obj.name,'train'))
                load('train','Fs','y');
            elseif (strcmp(obj.name,'splat'))
                load('splat','Fs','y');
            elseif (strcmp(obj.name,'laughter'))
                load('laughter','Fs','y'); 
            else
                cprintf('err', 'Error name must be chirp, gong, train, splat or laughter \n');return;
            end
            obj.signal = y;
        end

        function computeFourierRepresentation(obj)
              s.data = obj;
              fr = FourierTransformer.directTransform(s);
              obj.freq = fr;
        end

        function computeTimeRepresentation(obj)
              s.data = obj;
              fr = FourierTransformer.inverseTransform(s);
              obj.signal = fr;
        end        


    end
end

