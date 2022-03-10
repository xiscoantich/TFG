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
        S
        U
        V
    end
     
    properties (Access = private)
        name
    end
    
    methods (Access = public)

        function obj = Data(cParams)
            obj.init(cParams);
        end

        function plotSignal(obj)
           plot(obj.signal,'DisplayName','Signal')
           hold on;
           title('Signal');
           if(isempty(obj.rec_f) == false)
               plot(obj.rec_f,'DisplayName','Rec Fourier');
               hold on;
           end
           if(isempty(obj.rec_w) == false)
               plot(obj.rec_w,'DisplayName','Rec Wavelets');
               hold on;
           end
           if(isempty(obj.rec_pca) == false)
               plot(obj.rec_pca,'DisplayName','Rec PCA');
               hold on;
           end
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
                    obj.loadAudioSignal();
                    obj.dim = size(obj.signal,2);
                    obj.type_ft = cParams.type_ft;
                    obj.computeFourierRepresentation();
                    obj.computeWaveletRepresentation();
                    obj.computePCARepresentation()
                    
                case 'FREQUENCY'
                    obj.freq = cParams.freq;
                    obj.dim = size(obj.freq,2);
                    obj.type_ft = cParams.type_ft;
                    obj.computeTimeRepresentationFromFreq();
                    obj.wave = cParams.wave;
                    obj.wave_info = cParams.wave_info; 
                    obj.computeTimeRepresentationFromWave(); 
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
                    w = 15;
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
            obj.wave_info = wt;
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
            obj.rec_f = ift;
            %Aqui tengo que hacer para que tambien tenga la informacion de
            %que tipo de fft quiere hacer que deberia ser igual que la
            %anterior
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            wt = WaveletTransformer();
            iwt = wt.inverseTransform(s);
            obj.rec_w = iwt;
        end
        
        function computeTimeRepresentationFromPCA(obj)
            s.data = obj;
            pca = PCATransformer();
            ipca = pca.inverseTransform(s);
            obj.rec_pca = ipca;
        end
    end
end

