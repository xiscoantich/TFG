classdef FourierTransformer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        Data
    end
    
    methods (Access = public)
        function obj = FourierTransformer(Data)
            obj.Data = Data;
        end 
        function Transform(obj)
            %Determine if it is 1D or 2D
            %Compute Fourier Transform
            if isempty(obj.Data.signal)
                cprintf('err', 'FourierTransform error \n');return;  
            elseif length(obj.Data.signal(1,:))==1 || length(obj.Data.signal(:,1))==1
                obj.Data.ft = FFTCT_matrix(obj.Data.signal);
            else
                obj.Data.ft = FFT2CT(obj.Data.signal);
            end
        end
        
        function InverseTransform(obj)
            %Determine if it is 1D or 2D
            %Compute Inverse Fourier Transform
            if isempty(obj.Data.keep)
                obj.Data.keep = 1;
            end
            
            for i=1:1:length(obj.Data.keep)
                obj.Data.ft_cut = cut(obj.Data.ft, obj.Data.keep(i)); %Esto no se porque no funciona
                if isempty(obj.Data.ft)
                    cprintf('err', 'InverseTransform error \n');return;
                elseif length(obj.Data.signal(1,:))==1 || length(obj.Data.signal(:,1))==1
                    obj.Data.rec(:,:,i) = IFFTCT(obj.Data.ft_cut);
                else
                    obj.Data.rec(:,:,i) = IFFT2CT(obj.Data.ft_cut); 
                end          
            end
        end
    end
    
    methods (Access = private)
        function Atlow = cut(Bt, keep)
            Btsort = sort(abs(Bt(:)));
            thresh = Btsort(floor((1-keep)*length(Btsort)));
            ind = abs(Bt)>thresh;       %Find small index;
            Atlow = Bt.*ind;            %Theshold small indices
        end
    end
end