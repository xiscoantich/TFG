classdef Data < handle

    properties (Access = public)
        signal
        rec_pca
        dim
        freq
        sfreq
        wave
        motherwave
        dt
        wave_info
        type
        S
        U
        V
        Fs
        
        %New Wavelet Transformer
        level
        par
        ent_par
    end
     
    properties (Access = private)
        name
        originalsize
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
        
        function plotWaveCoef2(obj,mode)     
            %   Plot wavelet image (2D) decomposition.
            %   A short and simple function for displaying wavelet image decomposition
            %   coefficients in 'tree' or 'square' mode
            %
            %   Required : MATLAB, Image Processing Toolbox, Wavelet Toolbox
            %
            %   plotwavelet2(C,S,level,wavelet,rv,mode)
            %
            %   Input:  C : wavelet coefficients (see wavedec2)
            %           S : corresponding bookkeeping matrix (see wavedec2)
            %           level : level decomposition
            %           wavelet : name of the wavelet
            %           rv : rescale value, typically the length of the colormap
            %                (see "Wavelets: Working with Images" documentation)
            %           mode : 'tree' or 'square'
            %
            %   Output:  none
            %
            %   Example:
            %
            %     % Load image
            %     load wbarb;
            %     % Define wavelet of your choice
            %     wavelet = 'haar';
            %     % Define wavelet decomposition level
            %     level = 2;
            %     % Compute multilevel 2D wavelet decomposition
            %     [C S] = wavedec2(X,level,wavelet);
            %     % Define colormap and set rescale value
            %     colormap(map); rv = length(map);
            %     % Plot wavelet decomposition using square mode
            %     plotwavelet2(C,S,level,wavelet,rv,'square');
            %     title(['Decomposition at level ',num2str(level)]);
            %
            %
            %   Benjamin Tremoulheac, benjamin.tremoulheac@univ-tlse3.fr, Apr 2010
            if (nargin < 2) || isempty(mode), mode = 'square'; end
            
            C = obj.wave;
            S = obj.wave_info.l;
            wavelet = obj.motherwave;
            level = obj.level;
            
            %Define colormap and set rescale value
            colormap(gray); rv = length(gray);
            
            A = cell(1,level); H = A; V = A; D = A;
            for k = 1:level
                A{k} = appcoef2(C,S,wavelet,k); % approx
                [H{k} V{k} D{k}] = detcoef2('a',C,S,k); % details
                
                A{k} = wcodemat(A{k},rv);
                H{k} = wcodemat(H{k},rv);
                V{k} = wcodemat(V{k},rv);
                D{k} = wcodemat(D{k},rv);
            end
            if strcmp(mode,'tree') 
                aff = 0;  
                for k = 1:level
                    subplot(level,4,aff+1); image(A{k});
                    title(['Approximation A',num2str(k)]);
                    subplot(level,4,aff+2); image(H{k});
                    title(['Horizontal Detail ',num2str(k)]);
                    subplot(level,4,aff+3); image(V{k});
                    title(['Vertical Detail ',num2str(k)]);
                    subplot(level,4,aff+4); image(D{k});
                    title(['Diagonal Detail ',num2str(k)]);
                    aff = aff + 4;
                end
            elseif strcmp(mode,'square')
                dec = cell(1,level);
                dec{level} = [A{level} H{level} ; V{level} D{level}];
                for k = level-1:-1:1
                    dec{k} = [imresize(dec{k+1},size(H{k})) H{k} ; V{k} D{k}];
                end
                image(dec{1});
            end 
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            
            switch cParams.domain
                
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
                    obj.type = cParams.type;
                    obj.wave_info.level = cParams.level;
                    obj.wave_info.par = cParams.par;
                    obj.wave_info.ent_par = cParams.ent_par;
                    obj.motherwave = cParams.motherwave;
                    
                    obj.computeFourierRepresentation();
                    obj.computeWaveletRepresentation();
                    obj.computePCARepresentation()
                    
                case 'FOURIER'
                        obj.freq = cParams.freq;
                        obj.dim = cParams.dim;
                        obj.type.ft = cParams.type;
                        obj.originalsize = cParams.originalsize;
                        obj.computeTimeRepresentationFromFreq();
                  
                case 'WAVELET'
                        obj.wave = cParams.wave;
                        obj.dim = cParams.dim;
                        obj.type.wt = cParams.type;
                        obj.wave_info = cParams.wave_info;
                        obj.motherwave = cParams.motherwave;
                        obj.originalsize = cParams.originalsize;
                        %obj.level = cParams.level;
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
            switch obj.name
                case 'test'
                    obj.signal = magic(8);
                otherwise
                    path0=fileparts(which(mfilename));
                    imagepath=fullfile(path0,'Images',[obj.name,'.jpg']);
                    A = imread(imagepath);
                    obj.signal=double(rgb2gray(A));
            end
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
            obj.cutOriginalSize();
        end
        
        function computeTimeRepresentationFromWave(obj)
            s.data = obj;
            wt = WaveletTransformer();
            iwt = wt.inverseTransform(s);
            obj.signal = iwt;
            obj.cutOriginalSize();
        end
        
        function computeTimeRepresentationFromPCA(obj)
            s.data = obj;
            pca = PCATransformer();
            ipca = pca.inverseTransform(s);
            obj.signal = ipca;
        end
        
        function cutOriginalSize(obj)
            n = size(obj.signal);
            n0 = obj.originalsize;
            if mod(n0(1),2)~=0
                obj.signal = obj.signal(end-1,:);
                n0(1) = n0(1)-1;
            end
            if mod(n0(2),2)~=0
                obj.signal = obj.signal(:,end-1);
                n0(2) = n0(2)-1;
            end
            if n0(1) ~= n(1)
                ncut(1) = (n(1)-n0(1))/2;
                obj.signal = obj.signal(1+ncut(1):end-ncut(1),:);
            end
            if n0(2) ~= n(2)
                ncut(2) = (n(2)-n0(2))/2;
                obj.signal = obj.signal(:,1+ncut(2):end-ncut(2));
            end
        end
    end
end

