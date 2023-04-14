classdef Domain
    
    properties
        bound
    end
    
    methods (Abstract)        
        [x,dxdz] = reference2domain(obj, z)
        [z,dzdx] = domain2reference(obj, x)
        %
        [logdxdz,logdxdz2] = reference2domain_log_density(obj, z)
        [logdzdx,logdzdx2] = domain2reference_log_density(obj, x)
    end
    
    
    methods
        function y = domain(obj)
            y = obj.bound;
        end
        
        function y = domain_left(obj)
            y = obj.bound(1);
        end

        function y = domain_right(obj)
            y = obj.bound(2);
        end
    end
end