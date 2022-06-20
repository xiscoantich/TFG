classdef Data < handle
    properties (Access = public)
        signal
        audioinf
    end
    
    properties (Access = private)
        typesignal
        filename
        
    end
    methods (Access = public)
        function obj = Data(cParams)
            obj.init(cParams);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.typesignal = cParams.typesignal;
            obj.filename = cParams.filename;
            obj.loadData()
        end
        
        function obj = loadData(obj)
            switch obj.typesignal
                case {'audio','AUDIO', 'aud', 'AUD'}
                    obj.loadAudio();
                    
                case {'image','IMAGE', 'img', 'IMG'}
                    obj.loadImage();
                    
                otherwise
                    cprintf('err', 'typesignal err:  type audio or image  \n');
            end
        end
        
        function loadAudio(obj)
            switch obj.filename
                case {'chirp','gong','train','splat','laughter'}
                    load(obj.filename,'Fs','y');
                    obj.signal = y;
                    obj.audioinf.dt = 1/Fs;
                    obj.audioinf.Fs = Fs;
                case {'sinus'}
                    obj.audioinf.dt = 0.05;
                    t = 0:obj.audioinf.dt:1-obj.audioinf.dt;
                    w = 15;
                    obj.signal(:,1) = sin(w*t);
                
            end
            %Check if signal is a column vector.
            IsColumn = iscolumn(obj.signal);
            if IsColumn == 0 %False
                obj.signal = obj.signal.';
            end
         
        end
        
        function loadImage(obj)
            switch obj.filename
                case 'test'
                    obj.signal = magic(8);
                otherwise
                    path0=fileparts(which(mfilename));
                    imagepath=fullfile(path0,'/Images',[obj.filename]);
                    A = imread(imagepath);
                    obj.signal=double(rgb2gray(A));
            end
        end
            
    end
end