classdef Chebyshev2ndCDF < Chebyshev2nd & SpectralCDF
    
    methods
        function obj = Chebyshev2ndCDF(poly, varargin)
            obj@Chebyshev2nd(poly.order*2); 
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
            %bb = [0.5*sin(2*theta)-theta, sin((obj.n(2:end)+2).*theta)./(obj.n(2:end)+2) - sin(obj.n(2:end).*theta)./obj.n(2:end)]/pi;
            %
            cdf_ind = reshape(1:(obj.order+2), 1, []);
            tmp = sin(cdf_ind.*theta)./cdf_ind; % from 0 to n
            b = [tmp(:,2) - theta, tmp(:,3:end)-tmp(:,1:obj.order)]/pi;
            %
            %plonorm(b - bb, 'fro')
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            %
            theta = real(acos(x(:)));
            %
            cdf_ind = reshape(1:(obj.order+2), 1, []);
            tmp = sin(cdf_ind.*theta)./cdf_ind; % from 0 to n
            b = [tmp(:,2) - theta, tmp(:,3:end)-tmp(:,1:obj.order)]/pi;
            %
            % including the weight
            db = sin(theta.*(obj.n+1))*2/pi;
        end
        
    end
    
end