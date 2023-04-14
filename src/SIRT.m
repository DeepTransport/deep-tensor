classdef SIRT
    % SIRT class
    %
    % SIRT Properties:
    %   int_dir     - The direction for marginalising the approximation
    %                 >0: marginalise from x_d to x_1
    %                 <0: marginalise from x_1 to x_d
    %   fun_z       - Normalising constant of the function approximation part
    %   z           - Normalising constant
    %   tau         - denfensive term
    %   ref         - type of reference measure for the denfensive term
    %   oned_cdfs   - One dimensional bases for building CDFs.
    %
    % SIRT Methods:
    %   marginalise - Marginalise the squared FTT.
    %   eval_pdf    - Evaluate the normalised (marginal) pdf.
    %   eval_irt    - X = R^{-1}(Z), where X is the target random variable,
    %                 R is the Rosenblatt transport, and Z is the uniform
    %                 random variable.
    %               * This function can map marginal random variables.
    %   eval_cirt   - Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly
    %                 follow the target represented by SIRT, Z is uniform.
    %               * This function cannot handle marginal random variables.
    %   eval_rt     - Z = R(X), where Z is uniform and X is target.
    %               * This function can map marginal random variables.
    %   eval_rt_jac - Evaluate the Jacobian of Z = R(X).
    %               * This function cannot handle marginal random variables.
    %   random      - Generate random variables
    %   sobol       - Generate transformed sobol points
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also TTSIRT and PolySIRT
    
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
    end
    
    methods (Abstract)
        approximate(obj, func, base, opt, samples, deb)
        eval_irt_reference(obj, z);
        eval_rt_reference(obj, z);
        eval_potential_reference(obj, z);
        eval_cirt_reference(obj, z);
        eval_rt_jac_reference(obj, z);
        marginalise(obj, dir)
    end
    
    
    methods (Static)
        function y = get_potential_to_sqrt_density(base, y, z)
            [~,dxdz] = reference2domain(base, z);
            % log density of the reference measure
            mlogw = eval_measure_potential_reference(base, z);
            %
            % y is the potential function, use change of coordinates
            % y = exp( - 0.5*y + 0.5*sum(log(base.dxdz)) + 0.5*mlogw );
            y = exp( - 0.5*(y-mlogw-sum(log(dxdz),1)) );
        end
        
        function y = potential_to_sqrt_density(base, potential_fun, z)
            [x, dxdz] = reference2domain(base, z);
            y = feval(potential_fun, x);
            % log density of the reference measure
            mlogw = eval_measure_potential_reference(base, z);
            %
            % y is the potential function, use change of coordinates
            % y = exp( - 0.5*y + 0.5*sum(log(base.dxdz)) + 0.5*mlogw );
            y = exp( - 0.5*(y-mlogw-sum(log(dxdz),1)) );
        end
    end
    
    methods
        function obj = SIRT(defaultOpt, potential_func, base, varargin)
            defaultTau = 1E-8;
            defaultDebugger = Debugger();
            defaultSamples = [];
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            p = inputParser;
            %
            addRequired(p,'defaultOpt');
            addRequired(p,'potential_func');
            addRequired(p,'base',@(x) isa(x, 'ApproxBases') || isa(x, 'SIRT'));
            addOptional(p,'option',defaultOpt);
            %
            addParameter(p,'defensive', defaultTau, validErrTol);
            addParameter(p,'samples', defaultSamples);
            addParameter(p,'debug', defaultDebugger);
            addParameter(p,'previous', []);
            %
            p.KeepUnmatched = true;
            parse(p,defaultOpt,potential_func,base,varargin{:});
            opt = p.Results.option;
            tau = p.Results.defensive;
            deb = p.Results.debug;
            samples = p.Results.samples;
            previous = p.Results.previous;
            %
            if isa(base, 'SIRT')
                base = base.approx;
            end
            if ~isempty(samples)
                samples = domain2reference(base, samples);
            end
            if ~isempty(deb) && isevaluated(deb)
                z = domain2reference(base, deb.samples);
                deb = set_functions(deb, SIRT.get_potential_to_sqrt_density(base, deb.f, z));
                deb = set_samples(deb, z);
            end
            func = @(z) SIRT.potential_to_sqrt_density(base, potential_func, z);
            %
            obj.approx = approximate(obj, func, base, opt, samples, deb, previous);
            %
            obj.oned_cdfs = cell(size(obj.approx.oneds));
            for i = 1:length(obj.approx.oneds)
                obj.oned_cdfs{i} = CDFconstructor(obj.approx.oneds{i}, obj.approx.opt.cdf_tol);
            end
            obj.tau = tau;
            obj = marginalise(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = clean(obj)
            obj.approx = clean(obj.approx);
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
                [x,dxdr] = reference2domain(obj.approx, r, ind);
                %exp(-f) = exp(-f)/prod(obj.dxdz(ind));
                f = f + sum(log(dxdr),1);
            else
                [r,f,g] = eval_irt_reference(obj, z);
                [x,~] = reference2domain(obj.approx, r);
                [logdrdx,logdrdx2] = domain2reference_log_density(obj.approx, x);
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
            r = domain2reference(obj.approx, x, ind);
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
            bx = domain2reference(obj.approx, x, indx);
            [br,f] = eval_cirt_reference(obj, bx, z);
            [xr,dxdb] = reference2domain(obj.approx, br, indr);
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
            [r,drdx] = domain2reference(obj.approx, x);
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
            [r,drdx] = domain2reference(obj.approx, x, ind);
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