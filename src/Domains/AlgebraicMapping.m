classdef AlgebraicMapping < MappedDomain
    % forward mapping: from infinite (-\infty, \infty) to finite [-1,1]
    
    methods        
        function [z, dzdx] = domain2reference(obj, x)
            %[z, dzdx] = forward(obj, x)
            x = x/obj.scale;
            t = 1+x.^2;
            z = x.*t.^(-0.5);
            if nargout > 1
                dzdx = t.^(-1.5)/obj.scale;
            end
        end
        
        function [logdzdx,logdzdx2] = domain2reference_log_density(obj, x)
            x = x/obj.scale;
            t = 1+x.^2;
            logdzdx = -1.5*log(t)-log(obj.scale);
            if nargout > 1
                logdzdx2 = (-3/obj.scale)*(x./t);
            end
        end
        
        function [x, dxdz] = reference2domain(obj, z)
            %[x, dxdz] = inverse(obj, z)
            z(z<eps-1) = eps-1;
            z(z>1-eps) = 1-eps;
            t = 1-z.^2;
            t(t<eps) = eps;
            x = z.*t.^(-0.5)*obj.scale;
            if nargout > 1
                dxdz = t.^(-1.5)*obj.scale;
            end
        end
        
        function [logdxdz,logdxdz2] = reference2domain_log_density(obj, z)
            t = 1-z.^2;
            t(t<eps) = eps;
            logdxdz = -1.5*log(t) + log(obj.scale);
            if nargout > 1
                logdxdz2 = 3*(z./t);
            end
        end
    end
end