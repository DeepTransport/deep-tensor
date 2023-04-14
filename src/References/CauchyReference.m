classdef CauchyReference < SymmetricReference
    
    methods
        
        function  z = invert_ref_cdf(obj, u)
            z = tan((u-0.5)*pi);
        end
        
        function [u,f] = eval_ref_cdf(obj, z)
            u = atan(z)/pi + 0.5;
            f = 1./(pi*(1+z.^2));
        end
        
        function [f, g] = eval_ref_pdf(obj, z)
            f = 1./(pi*(1+z.^2));
            g = -(2*z)./(pi*(1+z.^2).^2);
        end
        
        function [f, g] = log_joint_ref_pdf(obj, z)
            f = - sum( log(1+z.^2), 1 ) - log(pi)*size(z,1);
            g = - (2*z)./(1+z.^2);
        end
        
        function obj = CauchyReference(varargin)
            obj@SymmetricReference(varargin{:});
        end
    end
    
end