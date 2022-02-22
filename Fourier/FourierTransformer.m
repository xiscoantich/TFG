classdef FourierTransformer < handle

    methods (Access = public, Static)

        function freq = directTransform(cParams)
            data = cParams.data;
            signal = data.signal;
            if data.dim == 1
                freq = FFTCT_matrix(signal);
            else
                freq = FFT2CT(signal);
            end
        end
        
        function signal = inverseTransform(cParams)
                data = cParams.data;
                freq = data.freq;
                if data.dim == 1
                    signal = IFFTCT(freq);
                else
                    signal = IFFT2CT(freq); 
                end
        end
    end


end