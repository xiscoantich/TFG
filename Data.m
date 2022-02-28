classdef Data < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        signal
        dim
        freq
        
        %Wavelet
        wave
        period
        scale
        coi
        dj
        paramout
        k
        
        %Esto hace falta?
        pad
        dj
        s0
        J1
        mother
        param
        
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
                    obj.computeTimeRepresentationFromFreq();
                    
                case 'TEMPORAL WAVE'
                    obj.name = cParams.name;
                    obj.loadAudioSignal();
                    obj.signal = cParams.signal;
                    obj.dt = cParams.dt;
                    obj.pad = cParams.pad;
                    obj.dj = cParams.dj;
                    obj.s0 = cParams.dj;
                    obj.J1 = cParams.J1;
                    obj.mother = cParams.mother;
                    obj.param = cParams.param;
                    obj.dim = size(obj.signal,2);
                    obj.computeWaveletRepresentation();
                    
                case 'FREQUENCY WAVE'
                    obj.wave = cParams.wave;
                    obj.dim = size(obj.freq,2);
                    obj.computeTimeRepresentationFromWave();
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
            ft = FourierTransformer.directTransform(s);
            obj.freq = ft;
        end
        
        function computeWaveletRepresentation(obj)
            s.data = obj;
            [wave,period,scale,coi, dj, paramout, k] = WaveletTransformer.directTransform(s);
            obj.wave = wave;
            obj.period = period;
            obj.scale = scale;
            obj.coi = coi;
            obj.dj = dj;
            obj.paramout = paramout;
            obj.k = k;
        end
        
        function computeTimeRepresentationFromFreq(obj)
            s.data = obj;
            ift = FourierTransformer.inverseTransform(s);
            obj.signal = ift;
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            iwt = WaveletTransformer.inverseTransform(s);
            obj.signal = iwt;
        end
        
    end
end

