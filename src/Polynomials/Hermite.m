classdef Hermite < Recurr

    methods
        function obj = Hermite(order)
            k = double((0:order))';
            a = ones(size(k));
            b = zeros(size(k));
            c = k;
            normalising = reshape( sqrt(1./cumprod(double([1, 1:order]))), 1, []);
            obj = obj@Recurr(order, a, b, c, normalising);
            obj.domain = [-inf,inf];
            obj.constant_weight = false;
        end
        
        function x = sample_measure(obj, n)
            x = randn(1,n);
        end        
        
        function x = sample_measure_skip(obj, n)
            x = sample_measure(obj, n);
        end
        
        function w = eval_measure(obj, x)
            w = exp(-0.5*x.^2)/sqrt(2*pi);
        end    

        function w = eval_log_measure(obj, x)
            w = -0.5*x.^2 - 0.5*log(2*pi);
        end   

        function w = eval_measure_deri(obj, x)
            w = -x.*exp(-0.5*x.^2)/sqrt(2*pi);
        end       
        
        function w = eval_log_measure_deri(obj, x)
            w = -x;
        end     
    end
    
end