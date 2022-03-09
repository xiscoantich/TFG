classdef Data < handle

    properties (Access = public)
        signal
        dim
        freq
        sfreq
        wave
        motherwave
        dt
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
        
        function plotWave(obj)
            imagesc(abs(obj.wave)) 
            xlabel('Time (integer index)') 
            ylabel('Scale')
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
                    obj.computeWaveletRepresentation();                    
                    
                case 'FREQUENCY'
                    obj.freq = cParams.freq;
                    obj.dim = size(obj.freq,2);
                    obj.computeTimeRepresentationFromFreq();
                    obj.wave = cParams.wave;
                    obj.computeTimeRepresentationFromWave();                    
                    
                    
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
                    load(obj.name,'Fs','y');
                    obj.signal = y;
                    obj.dt = 1/Fs;
                case {'sinus'}
                    obj.dt = 0.01;
                    t = 0:obj.dt:1;
                    w = 1000;
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
            wt = WaveletTransformer();
            [obj.wave] = wt.directTransform(s);
        end
        
        function computeTimeRepresentationFromFreq(obj)
            s.data = obj;
            ft = FourierTransformer();
            ift = ft.inverseTransform(s);
            obj.signal = ift;
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            wt = WaveletTransformer();
            iwt = wt.inverseTransform(s);
            obj.signal = iwt;
        end
        
    end
end

