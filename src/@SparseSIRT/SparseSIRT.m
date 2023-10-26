classdef SparseSIRT < SIRT 
    % SparseSIRT class
    %
    % SparseSIRT Properties:
    %   cum_P - The P matrix used for evaluating marginal density.
    %   batch_size  - 
    %           Batch size for multiple inputs
    %   bases_cdf_nodes - 
    %           Evaluation of basis for evaluating marginal density.
    %
    % SparseSIRT Methods:
    %   marginalise - 
    %           Marginalises the approximation.
    %   eval_pdf  - 
    %           Evaluates the normalised (marginal) pdf.
    %   eval_irt  - 
    %           X = R^{-1}(Z), where X is the target random variable, R is
    %           the Rosenblatt transport, and Z is the uniform random variable.
    %           * Cannot map marginal random variables.
    %   eval_cirt - 
    %           Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly follow
    %           the target represented by SIRT, Z is uniform.
    %           * This function cannot handle marginal random variables.
    %   eval_rt   - 
    %           Z = R(X), where Z is uniform and X is target.
    %           * Cannot map marginal random variables for sparse polynomial 
    %             approximations.
    %   eval_rt_jac - 
    %           Evaluates the Jacobian of Z = R(X).
    %           * This function cannot handle marginal random variables.
    %   random  - 
    %           Generates random variables
    %   sobol - Generates transformed sobol points
    %   set_defensive - 
    %           Resets the defensive term.
    %
    % See also SIRT, ONED, ONEDCDF, SparseOPTION, and SparseFun
    %
    
    properties
        bases_cdf_nodes
        cum_P
        batch_size
    end

    properties (Constant)
        defaultVar = InputData()
        defaultTau = 1E-8
        defaultOpt = SparseOption()
        defaultData = SparseData()
        defaultPoly = Legendre(30)
        defaultDomain = BoundedDomain([-1,1])
    end
    
    methods
        J = eval_rt_jac_reference_inner(obj, r, z)
    end

    methods        
        
        function obj = SparseSIRT(potential, base, varargin)
            obj@SIRT(potential, base, varargin{:});
            obj.batch_size = 1E2;
            obj.int_dir = 0;
        end

        function approx = approximate(obj, func, arg, opt, var)
            approx = SparseFun(func, arg, opt, 'var', var);
        end
    end
    
end