classdef LogarithmicMapping < MappedDomain
    % forward mapping: from infinite (-\infty, \infty) to finite [-1,1]
    
    methods
        function [z, dzdx] = domain2reference(obj, x)
            %[z, dzdx] = forward(obj, x)
            z = tanh(x/obj.scale);
            z(z<eps-1) = eps-1;
            z(z>1-eps) = 1-eps;
            if nargout > 1
                t = 1-z.^2;
                t(t<eps) = eps;
                dzdx = t/obj.scale;
            end
        end
        
        function [logdzdx,logdzdx2] = domain2reference_log_density(obj, x)
            z = tanh(x/obj.scale);
            z(z<eps-1) = eps-1;
            z(z>1-eps) = 1-eps;
            logdzdx = log(1-z.^2)-log(obj.scale);
            if nargout > 1
                logdzdx2 = (-2/obj.scale)*z;
            end
        end
        
        function [x, dxdz] = reference2domain(obj, z)
            %[x, dxdz] = inverse(obj, z)
            %z(z<eps-1) = eps-1;
            %z(z>1-eps) = 1-eps;
            x = real(atanh(z))*obj.scale;
            if nargout > 1
                t = 1-z.^2;
                t(t<eps) = eps;
                dxdz = obj.scale./t;
            end
        end
        
        function [logdxdz,logdxdz2] = reference2domain_log_density(obj, z)
            t = 1-z.^2;
            t(t<eps) = eps;
            logdxdz = -log(t)+log(obj.scale);
            if nargout > 1
                logdxdz2 = 2*(z./t);
            end
        end
        
    end
end