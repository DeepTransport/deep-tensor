classdef BoundedPolyCDF < Chebyshev2ndUnweighted & SpectralCDF
        
    methods
        function obj = BoundedPolyCDF(poly, varargin)
            obj@Chebyshev2ndUnweighted(poly.order*2); 
            obj@SpectralCDF(varargin{:});
        end

        function x = grid_measure(obj, n)
            x = linspace(obj.domain(1), obj.domain(2), n);
            x(1)   = obj.domain(1)-eps;
            x(end) = obj.domain(2)+eps;
        end

        function b = eval_int_basis(obj, x)
            %
            theta = real(acos(x));
            % the normalising constant used
            % int(U_n) = T_(n+1) / (n+1), n = 0, ..., order
            b = cos( theta.*(obj.n+1) ) .* (obj.normalising./(obj.n+1));
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            %
            % Evaluate Chebyshev polynomials of the 2nd kind,
            % for all input x, up to order n
            %
            % Inputs:
            % x:    n_pts
            %
            % Output:
            % f:    function outputs for each order at each x, n_pts x (n+1)
            % w:    weight function at each x, n_pts x 1
            %
            % x = cos(theta), or theta = acos(x), x in [-1, 1]
            %
            % U_n(theta) = cos( (n+1) * theta )  / cos(theta)
            % w(x) = sqrt(1 - x^2)
            %
            theta = real(acos(x(:)));
            b  = cos(theta.*(obj.n+1)) .* (obj.normalising./(obj.n+1));
            db = sin(theta.*(obj.n+1)) ./ (sin(theta)./obj.normalising);
            
            % deal with end points
            
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                db(mask,:) = repmat(((obj.n+1).*(-1).^obj.n).*obj.normalising, sum(mask), 1);
            end
            
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                db(mask,:) = repmat((obj.n+1).*obj.normalising, sum(mask), 1);
            end
        end
    end
    
end