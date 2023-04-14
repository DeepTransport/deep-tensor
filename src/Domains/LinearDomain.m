classdef LinearDomain < Domain
    
    properties
        mean
        dxdz
    end
    
    methods

        function y = shift(obj)
            y = obj.mean;
        end
        
        function y = scale(obj)
            y = obj.dxdz;
        end
        
        
        function [logdxdz,logdxdz2] = reference2domain_log_density(obj, z)
            logdxdz = log(obj.dxdz)*ones(size(z));
            logdxdz2 = zeros(size(z));
        end
        
        function [logdzdx,logdzdx2] = domain2reference_log_density(obj, x)
            logdzdx = (-log(obj.dxdz))*ones(size(x));
            logdzdx2 = zeros(size(x));
        end
        
        function [x,dxdz] = reference2domain(obj, z)
            x = z(:).*obj.dxdz+obj.mean;
            x = reshape(x,size(z));
            dxdz = obj.dxdz*ones(size(z));
        end
        
        function [z,dzdx] = domain2reference(obj, x)
            z = (x(:)-obj.mean)./obj.dxdz;
            z = reshape(z,size(x));
            dzdx = (1./obj.dxdz)*ones(size(x));
        end
    end
end