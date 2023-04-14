classdef ExpReference < HalfSpaceReference
    
    methods
        
        function  z = invert_ref_cdf(obj, u)
            tmp = 1-u;
            tmp(tmp < eps) = eps;
            z = -log(tmp);
        end
        
        function [u,f] = eval_ref_cdf(obj, z)
            u = 1-exp(-z);
            u(u<0) = 0;
            f = exp(-z);
        end
        
        function [f, g] = eval_ref_pdf(obj, z)
            f = exp(-z);
            g = -f;
        end
        
        function [f, g] = log_joint_ref_pdf(obj, z)
            f = -z;
            g = -ones(size(z));
        end
        
        function obj = ExpReference(varargin)
            obj@HalfSpaceReference(varargin{:});
        end
    end
    
end