classdef Chebyshev2ndTrigoCDF < Chebyshev2nd & TrigoCDF
    
    methods
        function obj = Chebyshev2ndTrigoCDF(poly, varargin)
            obj@Chebyshev2nd(poly.order*2); 
            obj@TrigoCDF(varargin{:});
        end

        function b = eval_int_basis(obj, theta)
            theta = theta(:);
            cdf_ind = reshape(1:(obj.order+2), 1, []);
            tmp = sin(cdf_ind.*theta)./cdf_ind; % from 0 to n
            b = [theta-tmp(:,2), tmp(:,1:obj.order)-tmp(:,3:end)]/pi;
        end
        
        function [b,db] = eval_int_basis_newton(obj, theta)
            theta = theta(:);
            cdf_ind = reshape(1:(obj.order+2), 1, []);
            tmp = sin(cdf_ind.*theta)./cdf_ind; % from 0 to n
            b = [theta-tmp(:,2), tmp(:,1:obj.order)-tmp(:,3:end)]/pi;
            %
            % including the weight
            db = (sin(theta.*(obj.n+1)).*sin(theta))*2/pi;
        end
        
    end
    
end