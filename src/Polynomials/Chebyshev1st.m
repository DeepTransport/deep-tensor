classdef Chebyshev1st < Spectral
    
    properties
        n(1,:) % only for Chebyshev
    end
    
    methods
        function obj = Chebyshev1st(order)
            obj.domain = [-1,1];
            obj.order = order;
            %
            n = obj.order + 1;
            obj.nodes = reshape( sort( double( cos( vpa(pi)*(2*(1:n)-1)/(2*n) ) ), 'ascend'), [], 1);
            obj.weights = double( ones(size(obj.nodes))/n );
            %
            obj.n = reshape(0:obj.order,1,[]);
            obj.normalising = reshape([1,sqrt(2)*ones(1,obj.order)], 1, []);
            %
            obj.constant_weight = false;
            obj = post_construction(obj);
        end
        
        function x = sample_measure(obj, n)
            z = rand(1,n); 
            x = sin(z*pi + pi/2);
        end   
        
        function w = eval_measure(obj, x)
            t = 1./(1-x.^2);
            t(t<eps) = eps;
            w = sqrt(t)/pi;
        end        

        function w = eval_log_measure(obj, x)
            t = 1 - x.^2;
            t(t<eps) = eps;
            w = -0.5*log(t) - log(pi);
        end   

        function w = eval_measure_deri(obj, x)
            t = 1-x.^2;
            t(t<eps) = eps;
            w = x.*t.^(-3/2)/pi;
        end       

        function w = eval_log_measure_deri(obj, x)
            t = 1-x.^2;
            t(t<eps) = eps;
            w = x./t;
            w(isinf(x)) = 1E16;
        end     
        
        function f = eval_basis(obj, x)
            %
            % Evaluate Chebyshev polynomials of the first kind,
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
            % T_n(x) = cos( n * theta )
            % w(x) = 1 / sqrt(1 - x^2)
            %
            theta = real(acos(x(:)));
            f = cos(theta.*obj.n) .* obj.normalising;
        end
        
        
        function f = eval_basis_deri(obj, x)
            %
            % Evaluate derivative of Chebyshev polynomials of the first kind,
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
            % T_n(x) = cos( n * theta )
            % w(x) = 1 / sqrt(1 - x^2)
            % d T_n(x) / dx = n * sin( n * theta) / sin(theta)
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
            if obj.order > 0
                f = [zeros(size(theta)), (sin(theta*obj.n(2:end)).*obj.n(2:end))./sin(theta)];
            else
                f = zeros(size(theta));
            end
            f = f.*obj.normalising;
        end
    end
end