classdef ApproxFun
    % ApproxFun class
    %
    % This class contains tensor-product polynomials basis functions and
    % mappings from the approximation domain to the reference domain.
    %
    % ApproxFun Properties:
    %   base  - See ApproxBases.
    %   data  - Data used for defining the approximation. 
    %           See also TTData and SparseData.
    %   opt   - Option that defines the approximation methods.
    %           See also TTOption and SparseOption.
    %   var   - User-specified sampling data for initializing TT and debugging
    %           various approximations. See also InputData.
    %   l2_err  - 
    %           l2 err computed by a cross-validation data set.
    %   l2_err  - 
    %           linf err computed by a cross-validation data set.
    %   n_eval  - 
    %           Number of function evaluations used.
    %
    % ApproxFun Methods:
    %   cardinals - 
    %           Cadinalities of the basis functions or the cardinality of
    %           the basis function for a chosen indexed coordinate.
    %   ndims - Dimension of the approximation domain.
    %   clean - Cleans up the approximation data and only keeps the results
    %           for evaluating the approximation.
    %   eval  - Evaluates the approximation. 
    %   grad  - Evaluates the gradient of the approximation. 
    %               
    % ApproxFun Methods (Abstract):
    %   approximate -
    %           Builds the approximation.
    %   eval_reference -
    %           Evaluates the approximation for reference variables.
    %   grad_reference -
    %           Evaluates the gradient for reference variables.
    %   int_reference  -
    %           Integrates the approximation for reference domain.
    %
    % See also TTFun and SparseFun.
    
    properties
        base ApproxBases
        data 
        opt  
        var
        l2_err
        linf_err
        n_eval
    end
    
    properties (Abstract, Constant)
        defaultVar
        defaultOpt
        defaultData
        defaultPoly
        defaultDomain
    end
    
    methods (Abstract)
        eval_reference(obj,z)
        grad_reference(obj,z)
        int_reference(obj)
    end
    
    methods (Static)
        function y = feval_reference(base, func, z)
            x = reference2domain(base, z);
            y = feval(func, x);
        end
    end
    
    methods
        
        function n = cardinals(obj,k)
            n = cardinals(obj.base,k);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function d = ndims(obj)
            d = ndims(obj.base);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = clean(obj)
            obj.data = clean(obj.data);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval(obj, x)
            % Evaluate the approximated function. f = EVAL(approx, x)
            %
            %   x   - input variables, d x n.
            %   f   - function values at x, 1 x n
            %
            
            z = domain2reference(obj.base, x);
            fx = eval_reference(obj, z);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [gx,fx] = grad(obj, x)
            % Evaluate the gradient of the approximation
            %   g = GRAD(approx, x)
            %
            %   x   - input variables, d x n, in the reference domain
            %   g   - gradient at x, d x n
            
            [z,dzdx] = domain2reference(obj.base, x);
            [gz, fx] = grad_reference(obj, z);
            gx = gz.*dzdx;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = isdebug(obj)
            flag = isdebug(obj.var);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [l2,linf] = rel_error(obj)
            if isdebug(obj.var)
                approx = eval_reference(obj, obj.var.debug_z);
                [l2,linf] = rel_error(obj.var, approx);
            else
                l2 = inf; 
                linf = inf;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = ApproxFun(func, arg, varargin)
            % Construct ApproxFun for a function mapping from R^d to R^m.
            %
            %   defaultOpt  - default ApproxFun options
            %   func - a function (R^d to R^m) that take inputs as a dxn
            %          matrix and returns mxn vector
            %   arg  - argument gives one of the following inputs
            %          (1) dimension of the input variable
            %          (2) a functional bases object
            %          (3) an existing ApproxFun used as the initial guess
            %          if (1) is chosen, then the fucntional bases is built
            %          using default domain and default basis
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parsing inputs
            
            p = inputParser;
            %
            addRequired (p,'func');
            addRequired (p,'arg');
            addOptional (p,'opt', obj.defaultOpt);
            addParameter(p,'var', obj.defaultVar);
            %
            p.KeepUnmatched = true;
            parse(p,func,arg,varargin{:});
            %
            obj.opt = p.Results.opt;
            obj.var = p.Results.var;
            %
            if isa(arg, 'ApproxFun')
                obj.base = arg.base;
                obj.data = arg.data; % init data
                obj.opt  = arg.opt;
            else
                if isa(arg, 'ApproxBases')
                    obj.base = arg;
                elseif isnumeric(arg) && isscalar(arg) && (arg > 0)
                    arg1 = obj.defaultPoly;
                    arg2 = obj.defaultDomain;
                    d = arg;
                    obj.base = ApproxBases(arg1,arg2,d);
                else
                    error('dimension should be a positive scalar')
                end
                obj.data = obj.defaultData; % init data
            end
            obj.var = set_debug(obj.var, func, obj.base);
            obj.n_eval = 0;
            obj.l2_err = inf;
            obj.linf_err = inf;
        end
    end
    
end