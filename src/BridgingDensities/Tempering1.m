classdef Tempering1 < SingleBeta
    
    methods
        function obj = Tempering1(varargin)
            defaultBetas = [];
            p = inputParser;
            addOptional(p,'betas',defaultBetas);
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            %
            obj@SingleBeta(varargin{:});
            obj.betas = p.Results.betas;
            %
            if isempty(obj.betas)
                obj.is_adapt = true;
                obj.betas = [];
            else
                obj.is_adapt = false;
                obj.betas = sort(obj.betas(:), 'ascend');
            end
            obj.n_layers = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_init(obj, mllkds, etol)
            if nargin < 3
                etol = 0.8;
            end
            etol = min(etol,obj.ess_tol_init);
            if obj.is_adapt
                % define beta_p
                beta_p = 0;
                beta = max(beta_p,obj.min_beta);
                % compute ess over sample size
                ess = ess_ratio((beta_p-beta)*mllkds);
                while ess > etol
                    beta = beta*obj.beta_factor;
                    ess = ess_ratio((beta_p-beta)*mllkds);
                end
                %out.beta = min(1, beta/obj.beta_factor);
                beta = min(1, beta);
                %
                log_weights = -beta*mllkds;
                log_weights = log_weights - max(log_weights);
                ess = ess_ratio(log_weights);
                %
                fprintf('\nInitial: beta=%3.3e, ess=%3.3e\n', beta, ess);
                %
                obj.init_beta = beta;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = adapt_density(obj, method, mllkds, mlps, mlogfx)
            if obj.is_adapt
                if obj.n_layers > 0
                    beta_p = obj.betas(obj.n_layers);
                    beta = max(beta_p,obj.min_beta);
                    switch method
                        case {'Aratio'}
                            % compute ess over sample size
                            ess = ess_ratio((beta_p-beta)*mllkds);
                            while ess > obj.ess_tol
                                beta = beta*obj.beta_factor;
                                ess = ess_ratio((beta_p-beta)*mllkds);
                            end
                        case {'Eratio'}
                            % compute ess over sample size
                            ess = ess_ratio(-beta*mllkds-mlps+mlogfx);
                            while ess > obj.ess_tol
                                beta = beta*obj.beta_factor;
                                ess = ess_ratio(-beta*mllkds-mlps+mlogfx);
                            end
                    end
                    %beta = min(1, beta/obj.beta_factor);
                    beta = min(1, beta);
                else
                    beta = obj.init_beta;
                end
                obj.betas = [obj.betas, beta];
            end
            obj.n_layers = obj.n_layers+1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [str, log_weights] = print_str(obj, mllkds, mlps, mlogfx)
            log_weights = -obj.betas(obj.n_layers)*mllkds-mlps+mlogfx;
            log_weights = log_weights - max(log_weights);
            ess = ess_ratio(log_weights);
            %
            % squared Hellinger error between logfx and (-beta_p*mllkds-mlps)
            if obj.n_layers > 1
                [~,dh2,~] = f_divergence(-mlogfx, -obj.betas(obj.n_layers-1)*mllkds-mlps);
                str = sprintf('\n\nIter=%2d, DHell=%3.3e, beta_p=%3.3e, beta=%3.3e, ess=%3.3e', ...
                    obj.n_layers, sqrt(dh2), obj.betas(obj.n_layers-1), obj.betas(obj.n_layers), ess);
            else
                str = sprintf('\n\nIter=%2d, beta=%3.3e, ess=%3.3e', ...
                    obj.n_layers, obj.betas(obj.n_layers), ess);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = islast(obj)
            flag = abs(obj.betas(obj.n_layers) - 1) < 1E-6;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = get_ratio_fun(obj, dirt, z, mllkds, mlps, ft) 
            beta = obj.betas(obj.n_layers);
            if obj.n_layers > 1 % not in the first layer
                % compute the reference density at z
                logfz = log_joint_pdf(dirt.ref, z);
                %
                beta_p = obj.betas(obj.n_layers-1);
                switch dirt.dirt_opt.method
                    case {'Aratio'}
                        f =  (beta-beta_p)*mllkds - logfz;
                    case {'Eratio'}
                        f =  beta*mllkds + mlps - ft - logfz;
                end
            else
                f = beta*mllkds + mlps;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = ratio_fun(obj, func, dirt, z, V) %(obj, func, beta, z, V)
            % Evaluate the ratio function during the DIRT construction
            %   f = RATIO_FUN(dirt, func, k, z, V)
            %
            %   func - function handle of the target function, returns minus log
            %          likelihood and minus log prior
            %   beta - the next beta value
            %   z    - reference random variables, d x n
            %   V    - dimension of the reduced basis or reparametrization basis
            %   f    - minus log pdf of the ratio function
            
            if nargin == 5 && ~isempty(V)
                [x, ft] = eval_irt(dirt, V*z);
            else
                [x, ft] = eval_irt(dirt, z);
            end
            % compute the minus log likelihood and minus log prior
            [mllkds, mlps] = func(x);
            f = get_ratio_fun(obj, dirt, z, mllkds, mlps, ft);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkds, mlps] = eval_wo_lis(obj, func, x)
            [mllkds, mlps] = func(x);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkds, mlps, mlogfx, gz_func, gz_ref, gz_irt] = ...
                eval_with_lis(obj, dirt, func, ratio_method, samples)
            if obj.n_layers == 0 % before adapt
                [mllkds, ~, gz_func] = func(samples);
                gz_irt = 0;
                gz_ref = 0;
                %
                mlps = zeros(1,size(samples,2));
                mlogfx = zeros(1,size(samples,2));
            else
                switch ratio_method
                    case {'Aratio'}
                        [x,mlogfx,~,Juz,Jux] = eval_irt(dirt, samples);
                        %
                        [mllkds, mlps, gx] = func(x);
                        gz_func = backtrack(dirt, Juz, Jux, gx);
                        gz_irt = 0;
                        gz_ref = 0;
                    case {'Eratio'}
                        [x,mlogfx,gz_irt,Juz,Jux] = eval_irt(dirt, samples);
                        %
                        [mllkds, mlps, gx] = func(x);
                        gz_func = backtrack(dirt, Juz, Jux, gx);
                        gz_ref = backtrack(dirt, Juz, Jux, x);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [V,s,g] = build_lis(obj, ratio_method, gz_func, gz_ref, gz_irt, log_weights)
            if obj.n_layers == 1 % post adapt
                beta = obj.betas(obj.n_layers);
                g = beta*gz_func;
            else
                beta_p = obj.betas(obj.n_layers-1);
                beta = obj.betas(obj.n_layers);
                switch ratio_method
                    case {'Aratio'}
                        g = (beta - beta_p)*gz_func;
                    case {'Eratio'}
                        g = beta*gz_func + gz_ref - gz_irt;
                end
            end
            weights = exp(log_weights)/sum(exp(log_weights));
            [V,s,~] = svd( g.*sqrt(weights(:)'), 'econ' );
            s = diag(s);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function y = restruct_ratio_fun(obj, dirt, arg, beta)
            
            if obj.n_layers > 1 % not in the first layer
                beta_p = obj.betas(obj.n_layers-1);
                if beta < beta_p
                    error('beta should be greater than beta_previous')
                end
                switch dirt.dirt_opt.method
                    case {'Aratio'}
                        rfval =  (beta-beta_p)*arg.mllkds - arg.logfz;
                    case {'Eratio'}
                        rfval =  beta*arg.mllkds + arg.mlps - arg.ft - arg.logfz;
                end
            else
                rfval = beta*arg.mllkds + arg.mlps;
            end
            
            y = exp( - 0.5*(rfval-arg.mlogw-arg.logJ) );
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = decompose_ratio_fun(obj, dirt, base, minus_log_prior, z, y, V)
            
            % y = exp( - 0.5*rfval + 0.5*sum(log(base.dxdz)) + 0.5*mlogw );
            mlogw = eval_measure_potential_reference(base, z);
            rfval = sum(log(base.dxdz)) + mlogw - log(y)*2;
            
            u = reference2domain(base, z);
            
            if nargin == 7 && ~isempty(V)
                [x, ft] = eval_irt(dirt, V*u);
            else
                [x, ft] = eval_irt(dirt, u);
            end
            
            % compute the reference density at z
            logfz = log_joint_pdf(obj.ref, u);
            
            % compute minus log prior values
            mlps = minus_log_prior(x);
            
            beta = obj.betas(obj.n_layers);
            if obj.n_layers > 1 % not in the first layer
                beta_p = obj.betas(obj.n_layers-1);
                switch dirt.dirt_opt.method
                    case {'Aratio'}
                        % rfval =  (beta-beta_p)*mllkds - logfz;
                        mllkds = (rfval + logfz)/(beta_p-beta);
                    case {'Eratio'}
                        % rfval =  beta*mllkds + mlps - ft - logfz;
                        mllkds = (rfval+ft+logfz-mlps)/beta;
                end
            else
                beta_p = -Inf;
                % rfval = beta*mllkds + mlps;
                mllkds = (rfval-mlps)/beta;
            end
            
            out.mllkds = mllkds;
            out.mlps = mlps;
            out.ft = ft;
            out.logfz = logfz;
            out.mlogw = mlogw;
            out.logJ = sum(log(base.dxdz));
            out.beta_p = beta_p;
            out.rfval = rfval;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function decompose_ratio_fun_debug(obj, dirt, base, minus_log_prior, func, z, y, V)
            %
            
            u = reference2domain(base, z);
            if nargin == 7 && ~isempty(V)
                x = eval_irt(dirt, V*u);
            else
                x = eval_irt(dirt, u);
            end
            
            beta = obj.betas(obj.n_layers);
            %
            out = decompose_ratio_fun(obj, dirt, base, minus_log_prior, z, y, V);
            
            % compute the minus log likelihood and minus log prior, for debug
            [mllkds_debug, mlps_debug] = func(x);
            
            disp(['log likelihood error: ' num2str(norm(mllkds_debug-out.mllkds)/norm(mllkds_debug))]);
            disp(['log prior error: ' num2str(norm(mlps_debug-out.mlps)/norm(mlps_debug))]);
            
            rf_debug = ratio_fun(obj, func, dirt, u, V);
            disp(['ratio fun error: ' num2str(norm(rf_debug-out.rfval)/norm(rf_debug))]);
            
            y_out = restruct_ratio_fun(obj, dirt, out, beta);
            disp(['func error: ' num2str(norm(y-y_out)/norm(y))]);
            
            beta2 = 0.5*(beta+out.beta_p);
            obj.betas(end) = beta2;
            rf_debug2 = ratio_fun(obj, func, dirt, u, V);
            y2 = exp( - 0.5*(rf_debug2-out.mlogw-out.logJ) );
            
            y_out2 = restruct_ratio_fun(obj, dirt, out, beta2);
            disp(['func error: ' num2str(norm(y2-y_out2)/norm(y2))]);
            
        end
        %}
        
    end
    
end