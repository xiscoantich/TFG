classdef Data < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        signal
        ft
        ft_cut
        rec
        fs
        keep
        id
        %fourier_transformer ???
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        function obj = Data()
            folder = fileparts(which(mfilename)); 
            % Add that folder plus all subfolders to the path.
            addpath(genpath(folder)); 
            
            %obj.fourier_transformer = FourierTransformer(obj); ???
        end
        
        function cut(obj)
           ft_sort =  sort(abs(obj.ft(:)));
           thresh = ft_sort(floor((1-obj.keep)*length(ft_sort)));
           ind = abs(obj.ft)>thresh;       %Find small index;
           obj.ft = obj.ft.*ind;            %Theshold small indices
        end
        
        function  imageread (obj)
           A = imread(obj.id);
           obj.signal=rgb2gray(A);
        end
        
        function audioread (obj)   
            [obj.signal, obj.Fs] = audioread(obj.id);
        end
        
        function audioload (obj)
            if (strcmp(obj.id,'chirp'))
                load('chirp','Fs','y'); obj.signal = y;
            elseif (strcmp(obj.id,'gong'))
                load('gong','Fs','y'); obj.signal = y;
            elseif (strcmp(obj.id,'train'))
                load('train','Fs','y'); obj.signal = y;
            elseif (strcmp(obj.id,'splat'))
                load('splat','Fs','y'); obj.signal = y;
            elseif (strcmp(obj.id,'laughter'))
                load('laughter','Fs','y'); obj.signal = y;
            else
                cprintf('err', 'Error id must be chirp, gong, train, splat or laughter \n');return;
            end
        end
    end
end

