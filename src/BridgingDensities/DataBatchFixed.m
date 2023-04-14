classdef DataBatchFixed < DataBatch
    
    methods
        function obj = DataBatchFixed(batches,varargin)
            p = inputParser;
            addRequired(p,'batches');
            p.KeepUnmatched = true;
            parse(p,batches,varargin{:});
            %
            obj@DataBatch(varargin{:});
            %
            batches = batches(:);
            obj.masks = cell(size(batches));
            obj.masks_inc = cell(size(batches));
            obj.data_order = [];
            for i = 1:length(batches)
                data_order_new = union(obj.data_order, batches{i}, 'stable');
                [~,inew] = setdiff(data_order_new, obj.data_order, 'stable');
                %
                obj.masks{i} = 1:numel(data_order_new);
                obj.masks_inc{i} = inew;
                %
                obj.data_order = data_order_new;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = set_init(obj, mllkds_vec, etol)
            if nargin < 3
                etol = 0.8;
            end
            etol = min(etol,obj.ess_tol_init);
            %
            ind = obj.data_order(obj.masks{1});
            mllkds = sum(mllkds_vec(ind,:),1);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = adapt_density(obj, method, mllkds_vec, mlps, mlogfx)
            if obj.data_layers == 0 && obj.n_layers == 0
                obj.data_layers = 1;
                %
                obj.betas = [0, obj.init_beta];
                obj.n_layers = obj.n_layers + 1;
                obj.beta_p = obj.init_beta;
            else
                if abs(obj.beta_p-1) < 1E-6
                    obj.data_layers = obj.data_layers + 1;
                    obj.beta_p = 0;
                end
                obj = adapt_temperature(obj, method, mllkds_vec, mlps, mlogfx);
            end
        end
        
    end
    
end