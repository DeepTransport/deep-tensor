classdef UniformReference < Reference
    
    properties
        pdf
    end
    
    methods
        
        function x = invert_cdf(obj, u)
            u(isinf(u)) = 1;
            u(isnan(u)) = 0;
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            x = (u./obj.pdf) + domain_left(obj.domain);
        end
        
        function [u,f] = eval_cdf(obj, x)
            u = (x - domain_left(obj.domain)) .* obj.pdf;
            f = ones(size(x))*obj.pdf;
        end
        
        function [f,g] = eval_pdf(obj, x)
            f = ones(size(x))*obj.pdf;
            g = zeros(size(x)); 
        end
        
        function [lf,glf] = log_joint_pdf(obj, x)
            lf = (log(obj.pdf)*size(x,1)) .* ones(1, size(x,2)); 
            glf = zeros(1, size(x,2)); 
        end
        
        function debug(obj, n)
            us = linspace(0, 1, n);
            xs = invert_cdf(obj, us);
            [c,f] = eval_cdf(obj, xs);
            fd = diff(c)/(10/n);
            figure;
            subplot(2,2,1)
            plot(xs(1:n-1),fd);hold on;plot(xs,f)
            legend('FD', 'pdf')
            title('pdf')
            subplot(2,2,2)
            plot(xs,c)
            title('cdf')
            %
            disp(sum(f)*(10/n))
            disp(norm(c-us));
            %
            [f,g] = eval_pdf(obj, xs);
            gd = diff(f)/(10/n);
            subplot(2,2,3)
            plot(xs(1:n-1),gd);hold on;plot(xs,g)
            legend('GD', 'grad')
            title('grad')
            subplot(2,2,4)
            plot(xs,f)
            title('pdf')
        end
        
        function z = random(obj, d, n)
            % pseudo random samples
            u = rand(d, n);
            z = invert_cdf(obj, u);
        end
        
        function z = sobol(obj, d, n)
            % QMC samples using Sobol sequence
            S = sobolset(d);
            u = net(S,n);
            z = invert_cdf(obj, u');
        end
        
        function obj = set_domain(obj, dom)
            obj.domain = dom;
            obj.pdf = 1/(domain_right(obj.domain) - domain_left(obj.domain));
        end
        
        function obj = UniformReference(dom)
            if (nargin<1)||(isempty(dom))
                dom = BoundedDomain([-1,1]);
            end
            obj = set_domain(obj, dom);
        end
        
    end
    
end