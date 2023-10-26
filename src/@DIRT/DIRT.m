classdef DIRT
    % DIRT class
    %
    % DIRT Properties:    
    %   ref         - The diagonal transformation and reference in DIRT
    %   logz        - Log normalising constant
    %   irts        - A cell array of IRTs
    %   method      - Type of DIRT construction, choose either 'Eratio'
    %                 or 'Aratio', default is 'Aratio'
    %   d           - dimension of the inpuit parameters
    %   n_layers    - Number of layers
    %   max_layers  - The maximum number of layers
    %   n_samples   - Number of samples used in the temperature adaptation 
    %               - and in the FTT construction
    %   n_debugs    - Number of debug samples
    %   adapt_beta  - A flag indicating if adaptive temperature is used
    %   min_beta    - Minimum temperature
    %   ess_tol     - Tolerance for increasing the temperature
    %   betas       - List of temperatures
    %   n_eval     - Number of function evaluations in TT cross.
    %
    % DIRT Methods:
    %   eval_irt    - X = R^{-1}(Z), where X is the target random variable,
    %                 R is the Rosenblatt transport, and Z is the uniform
    %                 random variable.
    %               * This function can map marginal random variables.
    %   eval_cirt   - Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly
    %                 follow the target represented by SIRT, Z is uniform.
    %               * This function cannot handle marginal random variables.
    %   eval_rt     - Z = R(X), where Z is uniform and X is target.
    %               * This function can map marginal random variables.
    %   backtrack   - Compute the gradient of f(T(z)) w.r.t. z
    %   ratio_fun   - The ratio function used in DIRT construction
    %   random      - Generate random variables
    %   sobol       - Generate transformed sobol points
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example: 
    %
    % % Use the following function in the example
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   function [mllkd, mlp] = fun_banana(u, data, sigma)
    %   F = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
    %   mllkd = sum((F-data).^2,1)/(2*sigma^2);
    %   mlp = 0.5*sum(u.^2,1);
    %   end
    %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Parameters for the target function:
    %   data=[3;5]; sigma=0.3;
    %
    % % Use the default option, Gaussian reference with 2nd order Langrange 
    % % basis, only need to provide the density funtion, dimension of
    % % random variables and domain
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,[-4,4]);
    %
    % % Case 1: using a uniform reference measure with Lagrange basis
    % %
    % % Define the map to the reference measure
    %   ref = UniformReference();
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Lagrange1(60,[-4,4]), Lagrange1(60,ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the approximate ratio method. The second argument
    % % is the dimenion of the parameters, 'ess_tol' is used for adaptation
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'min_beta',1E-3,'ess_tol',0.8,'method','Aratio');
    %
    % % Case 2: using a uniform reference measure with Legendre basis
    % %
    % % Define the map to the reference measure
    %   ref = UniformReference();
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Legendre(60, [-4, 4]), Legendre(40, ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the eaxct ratio method. Here the temperatures are
    % % prespecifed via the parameter 'betas'.
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'betas',2.^(-9:0),'method','Eratio');
    %
    % % Case 3: using a Gaussian reference measure with Fourier
    % %
    % % Define the map to the reference measure
    %   ref = GaussReference(0,1,[-5,5]);
    % % Define two basis functions, the first one is for discretising the 
    % % 0th level bridging density and the second one is for discretising 
    % % w.r.t. the reference measure.
    %   poly = {Fourier(30,[-5,5]), Fourier(30,ref.domain)};
    % % Define options for running ALS
    %   opt = FTToption('max_als',4,'init_rank',40,'max_rank', 50);
    % % Build DIRT using the approximate ratio method. The second argument
    % % is the dimenion of the parameters, 'ess_tol' is used for adaptation
    %   dirt = DIRT(@(u)fun_banana(u,data,sigma,beta),2,poly,ref,opt,...
    %       'min_beta',1E-3,'ess_tol',0.8,'method','Aratio');
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also SIRT
    
    properties
        logz
        irts
        %
        bridge
        %
        d
        ref
        dirt_opt
        %
        basis
        basis_r
        basis_s
        %
        lis_opt
        %
        n_eval
        %
        init_samples
        prev_approx
    end

    properties (Constant = true)
        defaultRef = GaussReference(0, 1, BoundedDomain([-4,4]));
        defaultDomain = BoundedDomain([-1,1])
        defaultPoly = Legendre(30)
        defaultDIRTOption = DIRTOption();
        defaultLISOption = LISOption();
        defaultInitSamples = [];
        defaultBridge = Tempering1('min_beta', 1E-3, 'ess_tol', 0.4);
    end
    
    properties (Abstract, Constant)
        defaultSIRTOpt
    end
    
    methods (Static)
        [bases,d] = build_bases(arg, ref)
        % build the approximation bases for DIRT
    end
    
    methods (Abstract)
        get_pre_sample_size(obj)
        get_inputdata(obj, base, n_layers, samples, density)
        get_new_layer(obj, func, bases, sirt_opt, n_layers, samples, density)
    end

    methods        
        [r,f,gz,Juz,Jux] = eval_irt(obj, z, k)
        % Evaluate DIRT r = T(z), where z follows some general reference
        
        [r,f] = eval_cirt(obj, x, z)
        % Evaluate the conditional DIRT
        
        [z,logz] = eval_rt(obj, r, k)
        % Evaluate deep RT z = T(r), where z is reference and r is target r.v.

        mlogz = eval_potential(obj, r, k)
               
        gz = backtrack(obj, Juz, Jux, gx)
        % Evaluate the gradient of f(T(z)), where T is a DIRT

        [f,g] = pullback(obj, func, z)
        % pullback of the target density function
        
        function n = num_layers(obj)    
            n = length(obj.irts);
        end
        
        function r = random(obj, n)
            % pseudo random samples
            z = random(obj.ref, obj.d, n);
            r = eval_irt(obj, z);
        end
        
        function r = sobol(obj, n)
            % QMC samples using Sobol sequence
            z = sobol(obj.ref, obj.d, n);
            r = eval_irt(obj, z);
        end

        function obj = DIRT(func, arg, varargin)
            % Call FTT constructor to build the FTT and setup data
            % structures for SIRT. Need to run marginalise after this.
            % parsing inputs
            %
            p = inputParser;
            addRequired(p, 'func',  @(x) isa(x, 'function_handle'));
            addRequired(p, 'arg');
            %
            addOptional(p, 'bridge', DIRT.defaultBridge);
            addOptional(p, 'ref', DIRT.defaultRef);
            %
            addOptional(p, 'sirt_option',obj.defaultSIRTOpt);
            addOptional(p, 'dirt_option',DIRT.defaultDIRTOption);
            addOptional(p, 'lis_option',DIRT.defaultLISOption);
            %
            addParameter(p, 'init_samples',DIRT.defaultInitSamples);
            addParameter(p, 'prev_approx',{});
            %
            p.KeepUnmatched = true;
            parse(p,func,arg,varargin{:});
            %
            obj.bridge = p.Results.bridge;
            obj.ref = p.Results.ref;
            sirt_opt = p.Results.sirt_option;
            obj.dirt_opt = p.Results.dirt_option;
            obj.lis_opt = p.Results.lis_option;
            obj.init_samples = p.Results.init_samples;
            obj.prev_approx = p.Results.prev_approx;
            %
            [bases, d] = DIRT.build_bases(arg, obj.ref);
            obj.d = d;
            %
            %
            switch obj.lis_opt.method
                case {'reduction'}    
                    obj.lis_opt = set_max_rank(obj.lis_opt, min(ceil(obj.lis_opt.rank_frac*d),d));
                    obj = build_lis(obj, func, bases, sirt_opt);
                case {'unitary'}    
                    obj = build_lis(obj, func, bases, sirt_opt);
                case {'none'}
                    obj = build(obj, func, bases, sirt_opt);
            end
        end
    end
    
end