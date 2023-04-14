classdef Reference
    % Reference class
    %               - Superclass for all one dimensional reference densities
    %
    % diagonalMap Methods:
    %   eval_cdf    - cdf function and pdf (2nd argument)
    %   eval_pdf    - pdf function and its gradient (2nd argument)
    %
    % See also UniformReference, GaussReference, LaplaceReference and
    % CauchyReference
    
    properties
        domain
    end
    
    methods (Abstract)
        [u,f] = eval_cdf(obj, x)
        [f,g] = eval_pdf(obj, x)
        %
        u = invert_cdf(obj, u)
        %
        [f,g] = log_joint_pdf(obj, x)
        %
        z = random(obj, d, n)
        z = sobol(obj, d, n)
        %
        obj = set_domain(obj, domain)
    end
    
    methods    
        function dom = get_domain(obj)
            dom = obj.domain;
        end
    end
    
end