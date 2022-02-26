classdef FourierTransformer < handle

    methods (Access = public, Static)

        function freq = directTransform(cParams)
            data = cParams.data;
            signal = data.signal;
            if data.dim == 1
                freq = FFT(signal);
            else
                freq = FFT2(signal);
            end 
        end
        
        function signal = inverseTransform(cParams)
                data = cParams.data;
                freq = data.freq;
                if data.dim == 1
                    signal = IFFT(freq);
                else
                    signal = IFFT2(freq); 
                end
        end
    end
    
    methods (Access = private)
        %Aqui deberia ir FFT, FFT2, IFFT, IFFT2
        %Sale un error y no se como solucionarlo
        %Como las funciones se utilizan unas a otras no pueden ser
        %estaticas, y sugen problemas 
    end
end