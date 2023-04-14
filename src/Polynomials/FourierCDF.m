classdef FourierCDF < Fourier & SpectralCDF
    
    methods
        function obj = FourierCDF(poly, varargin)
            obj@Fourier(poly.order*2); 
            obj@SpectralCDF(varargin{:});
        end

        function x = grid_measure(obj, n)
            x = linspace(obj.domain(1), obj.domain(2), n);
            x(1)   = obj.domain(1)-eps;
            x(end) = obj.domain(2)+eps;
        end

        function b = eval_int_basis(obj, x)
            %
            tmp = x(:).*obj.c;
            %b = [x(:), -cos(tmp)./(obj.c/sqrt(2)), sin(tmp)./(obj.c/sqrt(2))];
            b = [x(:), -cos(tmp)./(obj.c/sqrt(2)), sin(tmp)./(obj.c/sqrt(2)), ...
                sin(x(:)*(obj.m*pi))/(pi*obj.m/sqrt(2))];
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            %
            tmp = x(:).*obj.c;
            ct  = cos(tmp);
            st  = sin(tmp);
            %ct  = cos(tmp)/sqrt(0.5);
            %st  = sin(tmp)/sqrt(0.5);
            %
            %b = [x(:), -ct./(obj.c/sqrt(2)), st./(obj.c/sqrt(2))];
            b = [x(:), -ct./(obj.c/sqrt(2)), st./(obj.c/sqrt(2)), ...
                sin(x(:)*(obj.m*pi))/(pi*obj.m/sqrt(2))];
            %
            %db = [ones(size(x(:))), st*sqrt(2), ct*sqrt(2)];
            db = [ones(size(x(:))), st*sqrt(2), ct*sqrt(2), ...
                cos(x(:)*(obj.m*pi))*sqrt(2)];
        end
    end
    
end