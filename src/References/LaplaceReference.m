classdef LaplaceReference < SymmetricReference
    
    methods
        
        function  z = invert_ref_cdf(obj, u)
            ind1 = u < 0.5;
            ind2 = ~ind1;
            z = zeros(size(u));
            z(ind1) = log(u(ind1)*2);
            z(ind2) = -log((1-u(ind2))*2);
        end
        
        function [u,f] = eval_ref_cdf(obj, z)
            ind1 = z < 0;
            ind2 = ~ind1;
            u = zeros(size(z));
            u(ind1) = 0.5*exp(z(ind1));
            u(ind2) = 1 - 0.5*exp(-z(ind2));
            f = exp( -abs(z) )/2;
        end
        
        function [f, g] = eval_ref_pdf(obj, z)
            f = exp( -abs(z) )/2;
            g = -sign(z).*f;
        end
        
        function [f, g] = log_joint_ref_pdf(obj, z)
            f = - sum(abs(z),1) - log(2)*size(z,1);
            g = -sign(z);
        end
        
        function obj = LaplaceReference(varargin)
            obj@SymmetricReference(varargin{:});
        end
    end
    
end