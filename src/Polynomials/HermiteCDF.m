classdef HermiteCDF < Hermite & SpectralCDF
        
    methods
        function obj = HermiteCDF(poly, varargin)
            obj@Hermite(poly.order*2); 
            obj@SpectralCDF(varargin{:});
        end

        function x = grid_measure(obj, n)
            %{
            dom = cdf('Normal', [-6, 6], 0, 1);
            x = icdf('Normal', linspace(dom(1),dom(2),n), 0, 1);
            x = min(x, 6);
            x = max(x, -6);
            %}
            a = -10; %min(obj.nodes(1), -8);
            b = 10; %max(obj.nodes(end), 8);
            x = linspace(a,b,n);
        end

        function b = eval_int_basis(obj, x)
            b = eval_basis(obj, x);
            b(:,2:end) = -b(:,1:end-1).*(obj.normalising(2:end)./obj.normalising(1:end-1));
            w = eval_measure(obj, x(:));
            b = b.*w(:);
            b(:,1) = erf(x(:)/sqrt(2))/2; 
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            b = eval_basis(obj, x);
            w = eval_measure(obj, x(:));
            db = b.*w(:);
            %
            b(:,2:end) = -b(:,1:end-1)./obj.normalising(1:end-1).*obj.normalising(2:end);
            b = b.*w(:);
            b(:,1) = erf(x(:)/sqrt(2))/2; 
        end
    end
   

end