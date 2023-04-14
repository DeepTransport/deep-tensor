classdef Fourier < Spectral
    
    properties 
        m
        n
        is
        ic
        c
    end
    
    methods
        function obj = Fourier(order)
            obj.domain = [-1,1];
            obj.order = order;
            %
            n = obj.order*2 + 2;
            obj.nodes = reshape( sort( (2/n)*(1:n) - 1, 'ascend'), [], 1);
            obj.weights = ones(size(obj.nodes))/n;
            %
            obj.m   = obj.order+1;
            %obj.n   = obj.m*2;
            %obj.is  = 2:obj.m;
            %obj.ic  = (obj.m+1):(obj.n-1);
            obj.c   = reshape((1:obj.order)*pi, 1, []);
            %
            obj.constant_weight = true;
            obj = post_construction(obj);
            obj.node2basis(end,:) = obj.node2basis(end,:)/2;
        end
        
        
        function n = cardinal(obj)
            n = length(obj.nodes); %-1;
        end
        
        function x = sample_measure(obj, n)
            x = rand(1,n)*2 - 1;
        end        
        
        function w = eval_measure(obj, x)
            w = 0.5*ones(size(x));
        end        

        function w = eval_log_measure(obj, x)
            w = log(0.5)*ones(size(x));
        end
        
        function w = eval_measure_deri(obj, x)
            w = zeros(size(x));
        end       

        function w = eval_log_measure_deri(obj, x)
            w = zeros(size(x));
        end     
        
        function f = eval_basis(obj, x)
            %
            tmp = x(:).*obj.c;     
            %f = [ones(size(x(:))), sin(tmp)*sqrt(2), cos(tmp)*sqrt(2)];
            f = [ones(size(x(:))), sin(tmp)*sqrt(2), cos(tmp)*sqrt(2), ...
                cos( x(:)*(obj.m*pi) )*sqrt(2)];
        end
        
        function f = eval_basis_deri(obj, x)
            %
            tmp = x(:).*obj.c;
            %f = [zeros(length(x),1), cos(tmp).*(obj.c*sqrt(2)), -sin(tmp).*(obj.c*sqrt(2))];
            f = [zeros(length(x),1), cos(tmp).*(obj.c*sqrt(2)), -sin(tmp).*(obj.c*sqrt(2)), ...
                -sin(x(:)*(obj.m*pi))*(sqrt(2)*obj.m*pi)];
        end
    end
end