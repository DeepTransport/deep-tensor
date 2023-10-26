classdef Chebyshev2ndUnweighted < Spectral
    
    properties
        n(1,:) % only for Chebyshev
    end
    
    methods
        function obj = Chebyshev2ndUnweighted(order)
            obj.domain = [-1,1];
            obj.order = order;
            %
            n = obj.order + 1;
            obj.nodes = reshape( sort( double( cos( vpa(pi)*(1:n)/(n+1) ) ), 'ascend'), [], 1);
            obj.weights = double( sin( vpa(pi)*(1:n)/(n+1)).^2*2/(n+1) );
            %
            obj.n = reshape(0:obj.order,1,[]);
            obj.normalising = 1; %sqrt(2/pi); 
            obj.constant_weight = false;
            %
            obj = post_construction(obj);
        end
        
        function x = sample_measure(obj, n)
            x = betarnd(1.5,1.5,1,n);
            x = (2*x)-1;
        end   
        
        function x = sample_measure_skip(obj, n)
            left  = (min(obj.nodes) - 1)/2;
            right = (max(obj.nodes) + 1)/2;
            x = rand(1,n)*(right-left) + left;
        end

        function w = eval_measure(obj, x)
            t = 1 - x.^2;
            t(t<eps) = 0;
            w = 2*sqrt(t)/pi;
        end       
 
        function w = eval_log_measure(obj, x)
            t = 1 - x.^2;
            t(t<eps) = eps;
            w = 0.5*log(t) + log(2/pi);
        end   

        function w = eval_measure_deri(obj, x)
            t = 1./(1-x.^2);
            t(t<eps) = eps;
            w = -x.*sqrt(t)*2/pi;
        end       
        
        function w = eval_log_measure_deri(obj, x)
            error('not implemented')
        end     
        
        function f = eval_basis(obj, x)
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
            % deal with end points
            f = sin(theta.*(obj.n+1)) ./ (sin(theta)./obj.normalising);
            
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                f(mask,:) = repmat(((obj.n+1).*(-1).^obj.n).*obj.normalising, sum(mask), 1);
            end
            
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                f(mask,:) = repmat((obj.n+1).*obj.normalising, sum(mask), 1);
            end
        end
        
        function f = eval_basis_deri(obj, x)
            %
            % Evaluate derivative of Chebyshev polynomials of the second kind,
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
            % ( d theta / dx ) = inv(dx/dtheta) = inv( - sin(theta))
            %
            % sin(theta)^2 = (1-x^2)
            %
            % U_n(theta) = sin( (n+1) * theta )  / sin(theta)
            % w(x) = sqrt(1 - x^2)
            % d U_n(x) / dx = - (n+1) * cos((n+1)*theta) / sin(theta)^2 + sin((n+1)*theta) cos(theta) / sin(theta)^3
            %               = ( - (n+1) * cos((n+1)*theta) + x sin((n+1)*theta) / sin(theta) ) / (1-x^2)
            %               = ( (n+1) * cos((n+1)*theta) - x sin((n+1)*theta) / sin(theta) ) / (x^2-1)
            %
            theta = real(acos(x(:)));
            
            % deal with end points
            mask = abs(x+1) < eps;
            if sum(mask) > 0
                theta(mask) = pi;
            end
            mask = abs(x-1) < eps;
            if sum(mask) > 0
                theta(mask) = 0;
            end
            
            f = ( cos(theta*(obj.n+1)).*(obj.n+1) - sin(theta*(obj.n+1)).*(x(:)./sin(theta)) )./(x(:).^2-1);
            f = f.*obj.normalising;
        end
        
    end
    
end