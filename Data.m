classdef Data < handle

    properties (Access = public)
        signal
        dim
        freq
    end
    
    properties (Access = private)
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

        function plotSignal(obj)
           plot(obj.signal)
        end
    
        function plotFrequency(obj)
            plot(real(obj.freq))
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
                %    obj.computeWaveletRepresentation();                    
                    
                case 'FREQUENCY'
                    obj.freq = cParams.freq;
                    obj.dim = size(obj.freq,2);
                    obj.computeTimeRepresentationFromFreq();
              %      obj.computeWaveletRepresentation();                    
                    
                    
                case 'TEMPORAL WAVE'
%                     obj.name = cParams.name;
%                     obj.loadAudioSignal();
%                     obj.signal = cParams.signal;
%                     obj.dt = cParams.dt;
%                     obj.pad = cParams.pad;
%                     obj.dj = cParams.dj;
%                     obj.s0 = cParams.dj;
%                     obj.J1 = cParams.J1;
%                     obj.mother = cParams.mother;
%                     obj.param = cParams.param;
%                     obj.dim = size(obj.signal,2);
                    obj.wave = cParams.wave;
                    obj.dim = size(obj.wave,2);
                    obj.computeTimeRepresentationFromFreq();
                    obj.computeFourierRepresentation();
                    
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
            switch obj.name
                case {'chirp','gong','train','splat','laughter'}
                    obj.signal = y;
                case {'sinus'}
                    t = linspace(0,1,100);
                    w = 10;
                    obj.signal(:,1) = sin(w*t);
            end
        end

        function computeFourierRepresentation(obj)
            s.data = obj;
            ft = FourierTransformer();
            obj.freq = ft.directTransform(s);
        end
        
        function computeWaveletRepresentation(obj)
            s.data = obj;
            w = WaveletTransformer();
            [wave,period,scale,coi, dj, paramout, k] = w.directTransform(s);
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
            ft = FourierTransformer();
            ift = ft.inverseTransform(s);
            obj.signal = ift;
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            iwt = WaveletTransformer.inverseTransform(s);
            obj.signal = iwt;
        end
        
    end
end

