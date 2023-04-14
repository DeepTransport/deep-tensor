classdef HalfSpaceReference < Reference
    
    properties
        mu
        sigma
        left
        right
        is_truncated
    end
    
    methods (Abstract)
        [u,f] = eval_ref_cdf(obj, z)
        [f,g] = eval_ref_pdf(obj, z)
        z = invert_ref_cdf(obj, u)
        %
        [f,g] = log_joint_ref_pdf(obj, z)
    end
    
    methods
        
        function x = invert_cdf(obj, u)
            u(isinf(u)) = 1-eps;
            u(isnan(u)) = eps;
            u(u>(1-eps)) = 1-eps;
            u(u<eps) = eps;
            %
            u = u*(obj.right-obj.left) + obj.left;
            %
            z = invert_ref_cdf(obj, u);
            x = z*obj.sigma + obj.mu;
            if obj.is_truncated
                x = min(x, domain_right(obj.domain));
                x = max(x, domain_left (obj.domain));
            end
        end
        
        function [u,f] = eval_cdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [u,f] = eval_ref_cdf(obj, z);
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            %
            u = (u - obj.left)/(obj.right-obj.left);
            u = min(u, 1);
            u = max(u, 0);
            %
            f = f/obj.sigma/(obj.right-obj.left);
        end
        
        function[f,g] = eval_pdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [f,g] = eval_ref_pdf(obj, z);
            f = f/obj.sigma/(obj.right-obj.left);
            g = g/obj.sigma^2/(obj.right-obj.left);
        end
        
        function[f,g] = log_joint_pdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [f,g] = log_joint_ref_pdf(obj, z);
            f = f - log(obj.sigma)*size(x,1) - log((obj.right-obj.left))*size(x,1);
            g = g/obj.sigma;
        end
        
        function obj = set_domain(obj, dom)
            obj.domain = dom;
            if obj.is_truncated
                obj.left  = eval_ref_cdf(obj, (domain_left (dom)-obj.mu)/obj.sigma);
                obj.right = eval_ref_cdf(obj, (domain_right(dom)-obj.mu)/obj.sigma);
            else
                obj.left  = 0;
                obj.right = 1;
            end
        end
        
        function z = random(obj, d, n)
            % pseudo random samples
            u = rand(d, n);
            u = u*(obj.right-obj.left) + obj.left;
            z = invert_cdf(obj, u);
            if obj.is_truncated
                z(z >= domain_right(obj.domain)) = domain_right(obj.domain);
                z(z <= domain_left (obj.domain)) = domain_left (obj.domain);
            end
        end
        
        function z = sobol(obj, d, n)
            % QMC samples using Sobol sequence
            S = sobolset(d);
            u = net(S,n);
            u = u'*(obj.right-obj.left) + obj.left;
            z = invert_cdf(obj, u);
            if obj.is_truncated
                z(z >= domain_right(obj.domain)) = domain_right(obj.domain);
                z(z <= domain_left (obj.domain)) = domain_left (obj.domain);
            end
        end
        
        function debug(obj, n)
            xs = linspace(0, 10, n);
            [c,f] = eval_cdf(obj, xs);
            fd = diff(c)/(10/n);
            figure;
            subplot(3,2,1)
            plot(xs(1:n-1),fd);hold on;plot(xs,f)
            legend('FD', 'pdf')
            title('pdf')
            subplot(3,2,2)
            plot(xs,c)
            title('cdf')
            %
            disp(norm(invert_cdf(obj, c)-xs));
            %
            disp(sum(f)*(10/n))
            %
            [f,g] = eval_pdf(obj, xs);
            gd = diff(f)/(10/n);
            subplot(3,2,3)
            plot(xs(1:n-1),gd);hold on;plot(xs,g)
            legend('GD', 'grad')
            title('grad')
            subplot(3,2,4)
            plot(xs,f)
            title('pdf')
            
            [lf,lg] = log_joint_pdf(obj, xs);
            lgd = diff(lf)/(10/n);
            subplot(3,2,5)
            plot(xs(1:n-1),lgd);hold on;plot(xs,lg)
            legend('GD', 'grad')
            title('log grad')
            subplot(3,2,6)
            plot(xs,lf)
            title('log pdf')
        end
        

        function obj = HalfSpaceReference(varargin)
            defaultMu = 0;
            defaultSigma = 1;
            defaultDomain = BoundedDomain([0,10]);
            %
            p = inputParser;
            addOptional(p,'mu',defaultMu,@(x) isnumeric(x) && isscalar(x));
            addOptional(p,'sigma',defaultSigma,@(x) isnumeric(x) && isscalar(x) && (x>=0));
            addOptional(p,'domain',defaultDomain,@(x) isa(x,'Domain'));
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            %
            obj.mu = p.Results.mu;
            obj.sigma = p.Results.sigma;
            dom = p.Results.domain;
            obj.is_truncated = isa(dom, 'BoundedDomain');
            obj = set_domain(obj, dom);
        end
    end
    
end