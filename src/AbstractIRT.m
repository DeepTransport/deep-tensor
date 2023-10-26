classdef AbstractIRT
    % AbstractIRT class
    %
    % AbstractIRT Properties:
    %   int_dir - 
    %           The direction for marginalising the approximation
    %           >0: marginalise from x_d to x_1
    %           <0: marginalise from x_1 to x_d
    %   fun_z - Normalising constant of the function approximation part
    %   z     - Normalising constant
    %   tau   - denfensive term
    %   oned_cdfs - 
    %           One dimensional bases for building CDFs.
    %   approx  - 
    %           Approximation of the sqrt of the target density.
    %
    % AbstractIRT Methods:
    %   marginalise - 
    %           Marginalises the approximation.
    %   eval_pdf  - 
    %           Evaluates the normalised (marginal) pdf.
    %   eval_irt  - 
    %           X = R^{-1}(Z), where X is the target random variable, R is
    %           the Rosenblatt transport, and Z is the uniform random variable.
    %           * Can map marginal random variables for TT approximations.
    %           * Cannot map marginal random variables for sparse polynomial 
    %             approximations.
    %   eval_cirt - 
    %           Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly follow
    %           the target represented by SIRT, Z is uniform.
    %           * This function cannot handle marginal random variables.
    %   eval_rt   - 
    %           Z = R(X), where Z is uniform and X is target.
    %           * Can map marginal random variables for TT approximations.
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
    %%%%%%%%%%%%%%%%%
    %
    % see also SIRT and IRT
    
    properties
        int_dir
        order
        %
        fun_z
        z
        tau
        oned_cdfs
        %
        approx
        %
        approx_base
    end
    
    properties (Abstract, Constant)
        defaultVar
        defaultTau
        defaultOpt
        defaultData
        defaultPoly
        defaultDomain
    end
    
    methods (Abstract)
        get_potential_to_density(obj, base, y, z)
        potential_to_density(obj, base, potential_fun, z)
        %
        approximate(obj, func, base, opt, var)
        eval_irt_reference(obj, z);
        eval_rt_reference(obj, z);
        eval_potential_reference(obj, z);
        eval_cirt_reference(obj, z);
        eval_rt_jac_reference(obj, z);
        marginalise(obj, dir)
    end
    
    methods
        function [potential,arg,opt,var,tau,base] = parse_input(obj,potential,arg,varargin)
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            p = inputParser;
            %
            addRequired (p,'potential');
            addRequired (p,'arg',@(x) isnumeric(x) || isa(x, 'ApproxBases') ...
                || isa(x, 'ApproxFun') || isa(x, 'AbstractIRT'));
            addOptional (p,'opt',obj.defaultOpt);
            addParameter(p,'var',obj.defaultVar);
            addParameter(p,'defensive',obj.defaultTau,validErrTol);
            %
            p.KeepUnmatched = true;
            parse(p,potential,arg,varargin{:});
            opt = p.Results.opt;
            var = p.Results.var;
            tau = p.Results.defensive;
            %
            if isa(arg, 'AbstractIRT')
                base = arg.approx.base;
                arg  = arg.approx;
            elseif isa(arg, 'ApproxFun')
                base = arg.base;
            elseif isa(arg, 'ApproxBases')
                base = arg;
            else
                if isnumeric(arg) && isscalar(arg) && (arg > 0)
                    arg1 = obj.defaultPoly;
                    arg2 = obj.defaultDomain;
                    d = arg;
                    base = ApproxBases(arg1,arg2,d);
                    arg = base;
                else
                    error('dimension should be a positive scalar')
                end
            end
        end
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = clean(obj)
            obj.approx.data = clean(obj.approx.data);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [x,f,g] = eval_irt(obj, z)
            % Evaluate the inverse of the squared Rosenblatt transport 
            % X = R^{-1}(Z), where X is the target variable and Z is uniform.
            %   [X,f,g] = EVAL_IRT(irt, Z)
            %
            %   Z - uniform random variables, d x n
            %   X - random variable drawn form the pdf defined by SIRT
            %   f - pdf of X
            %   g - gradient of log pdf at X
            %
            d = ndims(obj.approx);
            dz = size(z,1);
            z = min(z, 1-1E-16);
            z = max(z, 1E-16);
            if nargout < 3
                [r,f] = eval_irt_reference(obj, z);
                if obj.int_dir > 0 % from left to right
                    ind = 1:dz;
                elseif obj.int_dir < 0 % from right to left
                    ind = ((d-dz)+1):d;
                elseif ~isempty(obj.order)
                    ind = obj.order(1:dz);
                else
                    error('either order or int_dir needs to be specified')
                end
                [x,dxdr] = reference2domain(obj.approx.base, r, ind);
                %exp(-f) = exp(-f)/prod(obj.dxdz(ind));
                f = f + sum(log(dxdr),1);
            else
                [r,f,g] = eval_irt_reference(obj, z);
                [x,~] = reference2domain(obj.approx.base, r);
                [logdrdx,logdrdx2] = domain2reference_log_density(obj.approx.base, x);
                f = f - sum(logdrdx,1);
                % need to revise g
                g = g.*exp(logdrdx) - logdrdx2;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_rt(obj, x)
            % Evaluate the squared Rosenblatt transport Z = R(X), where Z
            % is the uniform variable and X is the target variable.
            %   Z = EVAL_RT(irt, X)
            %
            %   Z - uniform random variables, d x n
            %   X - random variable drawn form the pdf defined by SIRT
            
            d = ndims(obj.approx);
            dr = size(x,1);
            if obj.int_dir > 0 % from left to right
                ind = 1:dr;
            elseif obj.int_dir < 0 % from right to left
                ind = ((d-dr)+1):d;
            elseif ~isempty(obj.order)
                ind = obj.order(1:dr);
            else
                error('either order or int_dir needs to be specified')
            end
            r = domain2reference(obj.approx.base, x, ind);
            z = eval_rt_reference(obj, r);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [xr,f] = eval_cirt(obj, x, z)
            % Evaluate the inverse of the conditional squared Rosenblatt
            % transport Y | X = R^{-1}(Z, X), where X is given, (X, Y)
            % jointly follow the target distribution represented by SIRT,
            % and Z is uniform.
            %   [Y,f] = EVAL_CIRT(irt, X, Z)
            %
            %   Z - uniform random variables, m x n
            %   X - given variables following the marginal of the target
            %   Y - random variable drawn from Y | X
            %   f - pdf of Y | X
            %
            % The conditioning depends on irt.int_dir.
            %   * >0 (marginalised from right to left), X = (X_1,...,X_k)
            %     is the the first k coordinates and Y = (X_{k+1},...,X_d)
            %   * <0 (marginalised from left to right), X = (X_m, ..., X_d)
            %     is the last (d-m+1) coordinates and Y = (X_1,...,X_{m-1})
            
            d = ndims(obj.approx);
            nr = size(z,2);
            dr = size(z,1);
            nx = size(x,2);
            dx = size(x,1);
            
            if dx == 0 || dr == 0
                error('dimension of x or the dimension of z should be nonzero')
            end
            
            if dr + dx ~= d
                error('dimension of x and the dimension of z mismatch the dimension of the joint')
            end
            
            if nx ~= nr && nx ~= 1
                error('number of x and the number of z mismatch')
            end
            
            % first compute the marginal
            if obj.int_dir > 0
                indx = 1:dx;
                indr = dx+(1:dr);
            elseif obj.int_dir < 0 
                indx = dr+(1:dx);
                indr = 1:dr;
            elseif ~isempty(obj.order)
                indx = obj.order(1:dx);
                indr = obj.order(dx+(1:dr));
            else
                error('either order or int_dir needs to be specified')
            end
            bx = domain2reference(obj.approx.base, x, indx);
            [br,f] = eval_cirt_reference(obj, bx, z);
            [xr,dxdb] = reference2domain(obj.approx.base, br, indr);
            %exp(-f) = exp(-f)/prod(obj.dxdz(indr));
            f = f + sum(log(dxdb),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function J = eval_rt_jac(obj, x, z)
            % Evaluate the jacobian of the squared Rosenblatt transport 
            % Z = R(X), where Z is the uniform random variable and X is the 
            % target random variable.
            %   J = EVAL_RT_JAC(irt, X, Z)
            %
            %   X - random variable drawn form the pdf defined by SIRT
            %   Z - uniform random variables, d x n
            %   J - Jacobian, d x (d x n), each d x d block is the Jabocian for X(:,j)
            
            % map to the reference coordinate
            [r,drdx] = domain2reference(obj.approx.base, x);
            J = eval_rt_jac_reference(obj, r, z); %dzdr
            % J - Jacobian, (d x d) x n, each d x d block is the Jabocian for X(:,j)
            % J = J.*dbdx(:);
            [d,n] = size(z);
            for i = 1:n
                ind = (i-1)*d + (1:d);
                J(:,ind) = J(:,ind).*drdx(:,i)';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval_pdf(obj, x)
            % Evaluate the normalised (marginal) pdf
            %   f = EVAL_PDF(irt, x)
            %
            %   x - input variables in the reference coordinates
            %   f - marginal density at x
            
            fx = eval_potential(obj, x);
            fx = exp(-fx);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval_potential(obj, x)
            % Evaluate the normalised (marginal) potential function
            %   f = EVAL_POTENTIAL(irt, x)
            %
            %   x - input variables in the reference coordinates
            %   f - marginal density at x
            
            % map to the reference coordinate
            
            d = ndims(obj.approx);
            dz = size(x,1);
            if obj.int_dir > 0
                ind = 1:dz;
            elseif obj.int_dir < 0 % from right to left
                ind = ((d-dz)+1):d;
            elseif ~isempty(obj.order)
                ind = obj.order(1:dz);
            end
            [r,drdx] = domain2reference(obj.approx.base, x, ind);
            fr = eval_potential_reference(obj, r);
            fx = fr - sum(log(drdx),1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = random(obj, n)
            % pseudo random samples
            d = length(obj.oned_cdfs);
            u = rand(d, n);
            r = eval_irt(obj, u);
        end
        
        function r = sobol(obj, n)
            % QMC samples using Sobol sequence
            d = length(obj.oned_cdfs);
            S = sobolset(d);
            u = net(S,n);
            r = eval_irt(obj, u');
        end
        
        function obj = set_defensive(obj, tau)
            obj.tau = tau;
            obj.z = obj.fun_z + tau;
        end
    end
    
end