classdef TTFun < ApproxFun
    % TTFun class
    %
    % This class builds and evaluates TT approximations to univariate and 
    % multivariate functions. 
    %
    %   - When building TT, it assumes that the function to be approximated 
    %     is defined on standarized intervals, e.g., [-1,1]^d for Chebyshev 
    %     polynomials.
    %
    %   - The TT approximation and its gradient can be evaluated in either 
    %     the reference coordinate or the orginal coordinates. The same 
    %     applies to integrating the TT.
    %
    %   - Outputs are horizontally aligned (dxn), where d is the dimension
    %     of the outputs and n is the number of input variables.
    %
    % TTFun Methods:
    %   cross       - Run the TT cross, this is called by the constructor.
    %
    %   The following do not consider the change of coordinate.
    %   eval_reference  - Evaluate TT. 
    %   grad_reference  - Evaluate the gradient of TT. 
    %   int_reference   - Integrate the entire TT.
    %
    %   The following consider the change of coordinate.
    %   eval        - Evaluate TT. Returns the same as eval_reference.
    %   grad        - Evaluate the gradient of TT. 
    %   int         - Integrate the entire TT.
    %
    %   Functions for internal usage.
    %   eval_block  - Evaluate TT for either the first or last k variables.
    %   int_block   - Integrate a block of TT cores.
    %
    %   Other functions.
    %   round       - Round the TT cores.
    %   size        - Size of the TT.
    %   rank        - Rank of the TT.
    
    %%%%%%%%%%%%%%%%%
    %
    % Example 1: (vector function outputs, m = 2):
    %
    % % Step 1: speficy the target function
    %   func = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
    %   d = 10; % dimensionality of the input
    %
    % % Step 2: setup the Legendre basis polynomial with order 20 in
    % % domain [0,1]
    %   poly = Legendre(20, [0,1]);
    %
    % % Step 3: use alternating energy enrichment (AMEN), default option
    %   opt = FTToption('max_als', 5, 'als_tol', 1E-8, 'local_tol', 1E-10, ...
    %           'kick_rank', 2, 'init_rank', 6, 'max_rank', 12);
    % % Optional: debug samples
    %   debug_size = 1E4;
    %   debug_x = zeros(d, debug_size);
    %   for k = 1:d
    %       debug_x(k,:) = sample_domain(poly, debug_size);
    %   end
    %
    % % Step 4: build FTT
    %   tt =  FTT(func, d, poly, opt, 'debug_x', debug_x);
    %
    % % Step 5: evaluate the function and its factorisation
    %   exact   = func(debug_x);
    %   appr_tt = eval(tt, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    % % Step 6: round the FTT by truncation local SVD. Here 1E-4 is the
    % % truncation threshold of each SVD (relative to the largest singular
    % % value)
    %   ttr = round(tt, 1E-4);
    %   appr_tt = eval(ttr, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example 2: (use the func and debug samples defined above)
    %
    % % Alternative Lagrange basis function, with order 5 and 4 elements
    %   poly = Lagrangep(5, 4, [0,1], 'ghost_size', 1E-10);
    %
    % % Alternative option use random enrichment
    %   opt_new = update(opt, 'tt_method', 'random');
    %
    % % Build FTT
    %   tt =  FTT(func, d, poly, opt_new, 'debug_x', debug_x);
    %
    % % Evaluate the function and its factorisation
    %   exact   = func(debug_x);
    %   appr_tt = eval(tt, debug_x);
    %   figure; plot(exact(:) - appr_tt(:), 'x')
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also ONED and FTTOPTION
    
    %{
    properties % from the super class
        data TTData
        base ApproxBases
        opt  TTOption
    end
    %}

    properties (Constant)
        defaultVar = InputData()
        defaultOpt = TTOption()
        defaultData = TTData()
        defaultPoly = Lagrangep(2,20)
        defaultDomain = BoundedDomain([-1,1])
    end
            
    methods (Static)
        [core, interp_x, res_w, res_x, core_next] = build_basis_amen(oned, ...
            interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr, ...
            dir, int_method, loc_err_tol, max_rank, kick_rank)
        % build tt core by amen
        
        [core,interp_x,core_next] = build_basis_svd(oned, ...
            interp_xold, core_next, F, ...
            dir, int_method, loc_err_tol, max_rank)
        % build tt core by svd
        
        [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F)
        % truncate local core
        
        interp_x = local_index(oned, direction, interp_xold, ind)
        % index selection
        
        rerr = local_error(core, f)
        % error of a core
        
        [f, f_evals] = local_block(oned, xleft, xright, func)
        % evaluate function for a coordinate
    end
    
    methods
        
        f = eval_block(obj, x, dir)
        
        obj = cross(obj, func, d, sample_x, debug_x)
        % cross iterations
        
        obj = round(obj, thres)
        % round the TT cores
        
        z = int(obj)
        % Integrate the entire TT
        
        % [g,f] = grad(obj, x)
        % gradient of the TT
        
        ftt = int_block(obj, ind)
        % Integrate a block of TT cores
        
        % fx = eval(obj, x, varargin)
        % Evaluate the ftt function. The output is horizontally aligned.
        
        % fx = eval_block(obj, x, dir)
        % Evaluate the fTT for either the first or last k variables.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function rs = rank(obj)
            rs = rank(obj.data);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = use_amen(obj)
            flag = strcmp(obj.opt.tt_method, 'amen');
        end
        
        function n = calc_sample_size(obj)
            n = obj.opt.init_rank + obj.opt.kick_rank*(1+obj.opt.max_als);
            n = n*ndims(obj.base);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = TTFun(func, arg, varargin)
            % Construct tensor train for a function mapping from R^d to R^m.
            %
            %   func - a function (R^d to R^m) that take inputs as a dxn
            %          matrix and returns mxn vector
            %   arg  - argument gives one of the following inputs
            %          (1) dimension of the input variable
            %          (2) a functional bases object
            %          (3) an existing TTFun used as the initial guess
            %          if (1) is chosen, then the fucntional bases is built
            %          using default domain and default basis
            %   opt  - TTFun options
            %
            obj@ApproxFun(func, arg, varargin{:});
            % add samples for initializing TT if needed
            nn = calc_sample_size(obj);
            obj.var = set_samples(obj.var, obj.base, nn);
            % update option: obj.opt.local_tol = obj.opt.local_tol/sqrt(ndims(obj)-1);
            if obj.opt.kick_rank == 0
                obj.opt.tt_method = 'fix_rank';
            end
            obj = cross(obj, func);
        end
    end
end