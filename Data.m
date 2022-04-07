classdef Data < handle

    properties (Access = public)
        signal
        rec_f
        rec_w
        rec_pca
        dim
        freq
        sfreq
        wave
        motherwave
        dt
        wave_info
        type_ft
        type_wt
        S
        U
        V
        Fs
        
        %New Wavelet Transformer
        level
    end
     
    properties (Access = private)
        name
    end
    
    methods (Access = public)

        function obj = Data(cParams)
            obj.init(cParams);
        end

        function plotSignal(obj,name)
           str = string(name);
           plot(obj.signal,'DisplayName',str)
           title('Signal');
           legend('show');
        end
    
        function plotFrequency(obj)
            plot(real(obj.freq))
            title('Frequency')
        end
        
        function plotWave(obj)
            imagesc(abs(obj.wave)) 
            xlabel('Time (integer index)') 
            ylabel('Scale')
        end
        
        function plotSurfWave(obj)
            surf(abs(obj.wave))
            hold on
            imagesc(abs(obj.wave))
        end
        
        function plotPCAInfo(obj)
            S = obj.S;
            subplot (1,2,1)
            semilogy(diag(S),'k','LineWidth',2)
            grid on
            xlabel('r')
            ylabel('Singular value,\sigma_r')
            %xlim([-50,1550])
            subplot(1,2,2)
            plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2)
            grid on
            xlabel('r')
            ylabel('Cumulative Energy')
            %xlim([-50,1550])
            ylim([0 1.1])
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            
            switch cParams.type
                
                case 'TEMPORAL'
                    obj.name = cParams.name;
                    switch cParams.typesignal
                        case 'AUDIO'
                            obj.loadAudioSignal();
                            obj.dim = 1;
                        case 'IMAGE'
                            obj.loadImage();
                            obj.dim = 2;
                        case 'VIDEO'
                            obj.loadVideo()
                            obj.dim = 3;
                    end
                    obj.type_ft = cParams.type_ft;
                    obj.type_wt = cParams.type_wt;
                    obj.level = cParams.level;
                    obj.motherwave = cParams.motherwave;
                    obj.computeFourierRepresentation();
                    obj.computeWaveletRepresentation();
                    obj.computePCARepresentation()
                    
                case 'FOURIER'
                        obj.freq = cParams.freq;
                        obj.dim = cParams.dim;
                        obj.type_ft = cParams.type_ft;
                        obj.computeTimeRepresentationFromFreq();
                  
                case 'WAVELET'
                        obj.wave = cParams.wave;
                        obj.dim = cParams.dim;
                        obj.type_wt = cParams.type_wt;
                        obj.wave_info = cParams.wave_info;
                        obj.motherwave = cParams.motherwave;
                        obj.computeTimeRepresentationFromWave();
             
                case 'PCA'
                        obj.dim = cParams.dim;
                        obj.U = cParams.U;
                        obj.S = cParams.S;
                        obj.V = cParams.V;
                        obj.computeTimeRepresentationFromPCA()
                    
                    
                    %Aqui estoy reescribiendo la reconstruccion de la señal
                    %Esto se deberia separar de en dos casos diferentes.
                    %El caso en el que reescribo a partir de la frecuencia
                    %y el caso en el que reescribo a partir de wavelet
                    
                    
                    
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
            path0=fileparts(which(mfilename));
            imagepath=[path0 '\Images\' obj.name '.jpg'];
            A = imread(imagepath);
            obj.signal=double(rgb2gray(A));
        end
        
        function loadAudioSignal(obj)
            switch obj.name
                case {'chirp','gong','train','splat','laughter'}
                    load(obj.name,'Fs','y');
                    obj.signal = y;
                    obj.dt = 1/Fs;
                    obj.Fs = Fs;
                case {'sinus'}
                    obj.dt = 0.05;
                    t = 0:obj.dt:1-obj.dt;
                    w = 15;
                    obj.signal(:,1) = sin(w*t);
                
            end
            %Check if signal is a column vector.
            IsColumn = iscolumn(obj.signal);
            if IsColumn == 0 %False
                obj.signal = obj.signal.';
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
            obj.wave_info.N = wt.N;
            obj.wave_info.scale = wt.scale;
            obj.wave_info.paramout = wt.paramout;
            obj.wave_info.k = wt.k;
            obj.wave_info.l = wt.l;
            obj.type_wt = wt.type_wt;
        end
        
        function computePCARepresentation(obj)
            s.data = obj;
            pca = PCATransformer();
            [obj.U,obj.S,obj.V] = pca.directTransform(s);
        end
        
        function computeTimeRepresentationFromFreq(obj)
            s.data = obj;
            ft = FourierTransformer();
            ift = ft.inverseTransform(s);
            obj.signal = ift;
            %Aqui tengo que hacer para que tambien tenga la informacion de
            %que tipo de fft quiere hacer que deberia ser igual que la
            %anterior
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            wt = WaveletTransformer();
            iwt = wt.inverseTransform(s);
            obj.signal = iwt;
        end
        
        function computeTimeRepresentationFromPCA(obj)
            s.data = obj;
            pca = PCATransformer();
            ipca = pca.inverseTransform(s);
            obj.signal = ipca;
        end
    end
end

