classdef Jacobi11 < Recurr
    
    methods
        function obj = Jacobi11(order)
            %
            k = double((0:order))';
            a = (2*k+3).*(k+2)./(k+1)./(k+3);
            b = zeros(size(k));
            c = (k+2)./(k+3);
            normalising = reshape( sqrt( (2*k+3).*(k+2)./(8*(k+1))*4/3 ), 1, []);
            obj@Recurr(order, a, b, c, normalising);
            obj.domain = [-1,1];
            obj.constant_weight = false;
            %
            %obj.weights = obj.weights*4/3;
            %obj.node2basis = obj.node2basis*4/3;
        end
        
        function x = sample_measure(obj, n)
            x = betarnd(2,2,1,n);
            x = (2*x)-1;
        end        
        
        function x = sample_measure_skip(obj, n)
            left  = (min(obj.nodes) - 1)/2;
            right = (max(obj.nodes) + 1)/2;
            x = rand(1,n)*(right-left) + left;
        end
        
        function w = eval_measure(obj, x)
            w = (1-x.^2)*3/4;
        end  

        function w = eval_log_measure(obj, x)
            w = log(1-x.^2) + log(3/4);
        end   

        function w = eval_measure_deri(obj, x)
            w = -x*3/2;
        end       
        
        function w = eval_log_measure_deri(obj, x)
            error('not implemented')
        end     
    end
end