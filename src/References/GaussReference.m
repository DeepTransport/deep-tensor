classdef GaussReference < SymmetricReference
    
    methods
        
        function  z = invert_ref_cdf(obj, u)
            z = erfinv(u*2 - 1)*sqrt(2);
        end
        
        function [u,f] = eval_ref_cdf(obj, z)
            u = 0.5*( 1 + erf(z/sqrt(2)) );
            f = exp( -0.5*z.^2 - 0.5*log(2*pi) );
        end
        
        function [f, g] = eval_ref_pdf(obj, z)
            f = exp( -0.5*z.^2 - 0.5*log(2*pi) );
            g = -z.*f; 
        end
        
        function [f, g] = log_joint_ref_pdf(obj, z)
            f = sum(- 0.5*z.^2,1) - 0.5*log(2*pi)*size(z,1);
            g = - z; 
        end
        
        function obj = GaussReference(varargin)
            obj@SymmetricReference(varargin{:});
        end
    end
    
end