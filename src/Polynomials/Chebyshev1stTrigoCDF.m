classdef Chebyshev1stTrigoCDF < Chebyshev1st & TrigoCDF

    methods
        function obj = Chebyshev1stTrigoCDF(poly, varargin)
            obj@Chebyshev1st(poly.order*2); 
            obj@TrigoCDF(varargin{:});
        end
        
        function b = eval_int_basis(obj, theta)
            theta = theta(:);
            if obj.order > 0
                b = [theta/pi, sin(obj.n(2:end).*theta).*((sqrt(2)/pi)./obj.n(2:end))];
            else
                b = theta/pi;
            end
        end
        
        function [b,db] = eval_int_basis_newton(obj, theta)
            theta = theta(:);
            if obj.order > 0
                b = [theta/pi, sin(obj.n(2:end).*theta).*((sqrt(2)/pi)./obj.n(2:end))];
            else
                b = theta/pi;
            end
            %
            db = cos(theta.*obj.n) .* obj.normalising / pi;
        end
    end
end