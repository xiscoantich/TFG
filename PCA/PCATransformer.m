classdef PCATransformer < handle
    properties (Access = public)
        S
        V
        U
    end
    
     methods (Access = public)
         function [U,S,V] = directTransform(obj,cParams)
            data = cParams.data;
            signal = data.signal;
            [U,S,V] = svd(signal,"econ");
         end
         
         function Xrec = inverseTransform(obj,cParams)
             data = cParams.data;
             U = data.U;
             S = data.S;
             V = data.V;
             Xrec = U*S*V';
         end
     end

end