classdef PCATransformer < handle
    properties (Access = public)
        transtype
        transmethod
        S
        V
        U
        originalsize
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        function obj = PCATransformer(cParams)
            obj.init(cParams);
        end
        
        function obj = directTransform(obj,data)
            [obj.U,obj.S,obj.V] = svd(data.signal,"econ");
        end
        
        function rec = inverseTransform(obj,data)
             U = data.U;
             S = data.S;
             V = data.V;
             rec = U*S*V';
        end
    end
    
     methods (Access = private, Static)
        
        
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.transtype = cParams.transtype;
            obj.transmethod = cParams.transmethod;
            %obj.originalsize = size(data.signal);
        end
    end
end