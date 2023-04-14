classdef LaguerreCDF < Laguerre & SpectralCDF
        
    methods
        function obj = LaguerreCDF(poly, varargin)
            obj@Laguerre(poly.order*2); 
            obj@SpectralCDF(varargin{:});
        end
        
        function x = grid_measure(obj, n)
            %{
            dom = cdf('exp', [1E-6, 20], 1);
            x = icdf('exp', linspace(dom(1),dom(2),n), 1);
            %}
            b = max(obj.nodes(end), 15);
            x = linspace(1E-16, b, n);
        end

        function b = eval_int_basis(obj, x)
            b = eval_basis(obj, x);
            w = eval_measure(obj, x(:));
            t = cumsum(b(:,1:end-1),2);
            b(:,2:end) = ( t.*(w(:).*x(:)) ) ./ (1:obj.order);
            b(:,1) = -w(:);
        end
        
        function [b,db] = eval_int_basis_newton(obj, x)
            b = eval_basis(obj, x);
            w = eval_measure(obj, x(:));
            db = b.*w(:);
            %
            t = cumsum(b(:,1:end-1),2);
            b(:,2:end) = ( t.*(w(:).*x(:)) ) ./ (1:obj.order);
            b(:,1) = -w(:);
        end
    end
    
end