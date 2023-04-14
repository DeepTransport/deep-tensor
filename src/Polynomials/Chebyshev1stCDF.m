classdef Chebyshev1stCDF < Chebyshev1st & SpectralCDF

    methods
        function obj = Chebyshev1stCDF(poly, varargin)
            obj@Chebyshev1st(poly.order*2); 
            obj@SpectralCDF(varargin{:});
        end
        
        function x = grid_measure(obj, n)
            x = linspace(obj.domain(1), obj.domain(2), n);
            x(1)   = obj.domain(1)-eps;
            x(end) = obj.domain(2)+eps;
        end

        function b = eval_int_basis(obj, x)
            %
            theta = real(acos(x(:)));
            if obj.order > 0
                b = [theta/pi, sin(obj.n(2:end).*theta).*((sqrt(2)/pi)./obj.n(2:end))];
            else
                b = theta/pi;
            end
            b = -b;
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            %
            theta = real(acos(x(:)));
            if obj.order > 0
                b = [theta/pi, sin(obj.n(2:end).*theta).*((sqrt(2)/pi)./obj.n(2:end))];
            else
                b = theta/pi;
            end
            b = -b;
            %
            db = cos(theta.*obj.n) .* obj.normalising;
            w = eval_measure(obj, x(:));
            db = db.*w;
        end
    end
end