function obj = cross(obj, func)
% Cross iterations for building the TT, only call this function if you know
% what you are doing

d = ndims(obj.base);

% determine dimension of output and prepare data
y = func(sample_measure_reference(obj.base,1));
obj.data.nout = numel(y);

% initilise TT cores
if isempty(obj.data.cores)
    % copy properties
    obj.data.direction = 1;
    % nested interpolation points
    obj.data.cores = cell(d,1);
    obj.data.interp_x = cell(d,1);
    %
    for k = d:-1:1
        tmp = get_samples(obj.var, obj.opt.init_rank);
        obj.data.interp_x{k} = tmp(k:d,:);
    end
    % interpolation weights
    obj.data.cores{1} = rand(1, cardinal(obj.base.oneds{1}), obj.opt.init_rank, obj.data.nout);
    for k = 2:d-1
        obj.data.cores{k} = rand(obj.opt.init_rank, cardinal(obj.base.oneds{k}), obj.opt.init_rank);
    end
    obj.data.cores{d} = rand(obj.opt.init_rank, cardinal(obj.base.oneds{d}), 1);
else
    obj.data.direction = -obj.data.direction;
end
% initialise the residual coordinates for AMEN
if use_amen(obj) && isempty(obj.data.res_x)
    obj.data.res_x = cell(d,1);
    % nested interpolation points for res
    if obj.data.direction > 0 % direction has already been flipped
        for k = d:-1:1
            tmp = get_samples(obj.var, obj.opt.kick_rank);
            obj.data.res_x{k} = tmp(k:d,:);
        end
    else
        for k = d:-1:1
            tmp = get_samples(obj.var, obj.opt.kick_rank);
            obj.data.res_x{k} = tmp(1:k,:);
        end
    end
end
% reinitialise the residual blocks for AMEN
if use_amen(obj) && isempty(obj.data.res_w)
    obj.data.res_w = cell(d,1);
    if obj.data.direction < 0 % direction has already been flipped
        for k = 1:(d-1)
            obj.data.res_w{k} = ones(obj.opt.kick_rank, size(obj.data.cores{k},3));
        end
        obj.data.res_w{d} = ones(size(obj.data.cores{d},1), obj.opt.kick_rank);
    else
        obj.data.res_w{1} = ones(obj.opt.kick_rank, size(obj.data.cores{1},3));
        for k = 2:d
            obj.data.res_w{k} = ones(size(obj.data.cores{k},1), obj.opt.kick_rank);
        end
    end
end

%
if isdebug(obj.var)
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  max_debug_E  mean_debug_E  cum#fevals\n');
else
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  cum#fevals\n');
end

% start
als_iter = 0;
rs = ones(1,d);
while true % run ALS
    if obj.data.direction > 0
        ind = 1:(d-1);
    else
        ind = d:-1:2;
    end
    %
    switch obj.opt.tt_method
        case {'fix_rank'}
            % start
            errs = zeros(1,d);
            for k = ind
                if obj.data.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.data.interp_x{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, Jx_left, obj.data.interp_x{k+1}, func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    obj.n_eval = obj.n_eval + nf;
                    % Schmidt operator and interpolation, k - 1 is the previous index
                    % push couple matrix to the right
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k+1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                        Jx_left, obj.data.cores{k+1}, F, obj.data.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k) = size(obj.data.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                    else
                        Jx_right = obj.data.interp_x{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, obj.data.interp_x{k-1}, Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    obj.n_eval = obj.n_eval + nf;
                    % Schmidt operator and interpolation, k + 1 is the previous index
                    % push couple matrix to the right
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k-1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                        Jx_right, obj.data.cores{k-1}, F, obj.data.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k-1) = size(obj.data.cores{k}, 1);
                end
            end
        case {'random'}
            % at the head, update the random enrichment set
            % enrich = get_samples(obj.var, obj.opt.kick_rank);
            % start
            errs = zeros(1,d);
            for k = ind
                if obj.data.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.data.interp_x{k-1};
                    end
                    enrich = get_samples(obj.var, obj.opt.kick_rank);
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, Jx_left, obj.data.interp_x{k+1}, func);
                    [Fe,ne] = TTFun.local_block(obj.base.oneds{k}, Jx_left, enrich(k+1:d,:),    func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    obj.n_eval = obj.n_eval + nf + ne;
                    % Schmidt operator and interpolation, k - 1 is the previous index
                    % push couple matrix to the right
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k+1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                        Jx_left, obj.data.cores{k+1}, cat(3,F,Fe), obj.data.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k) = size(obj.data.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                    else
                        Jx_right = obj.data.interp_x{k+1};
                    end
                    enrich = get_samples(obj.var, obj.opt.kick_rank);
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, obj.data.interp_x{k-1}, Jx_right, func);
                    % different from left iteration, k + 1 is the previous index
                    [Fe,ne] = TTFun.local_block(obj.base.oneds{k}, enrich(1:k-1,:),   Jx_right,  func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    obj.n_eval = obj.n_eval + nf + ne;
                    % Schmidt operator and interpolation, k + 1 is the previous index
                    % push couple matrix to the right
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k-1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                        Jx_right, obj.data.cores{k-1}, cat(1,F,Fe), obj.data.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k-1) = size(obj.data.cores{k}, 1);
                end
            end
        case {'amen'}
            % start
            errs = zeros(1,d);
            for k = ind
                if obj.data.direction > 0
                    if k == 1
                        Jx_left = [];
                        Jr_left = [];
                        Rw_left = 1;
                    else
                        Jx_left = obj.data.interp_x{k-1};
                        Jr_left = obj.data.res_x{k-1};
                        Rw_left = obj.data.res_w{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, Jx_left, obj.data.interp_x{k+1}, func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = TTFun.local_block(obj.base.oneds{k}, Jr_left, obj.data.res_x{k+1}, func);
                    obj.n_eval = obj.n_eval + nf + nr;
                    % evaluate update function at x_k nodes
                    if k > 1
                        [Fu,nu] = TTFun.local_block(obj.base.oneds{k}, Jx_left, obj.data.res_x{k+1}, func);
                        obj.n_eval = obj.n_eval + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.res_w{k}, obj.data.res_x{k}, obj.data.cores{k+1}] = TTFun.build_basis_amen(...
                        obj.base.oneds{k}, Jx_left, Jr_left, Rw_left, obj.data.res_w{k+1}, obj.data.cores{k+1}, F, Fu, Fr, ...
                        obj.data.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k) = size(obj.data.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                        Jr_right = [];
                        Rw_right = 1;
                    else
                        Jx_right = obj.data.interp_x{k+1};
                        Jr_right = obj.data.res_x{k+1};
                        Rw_right = obj.data.res_w{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = TTFun.local_block(obj.base.oneds{k}, obj.data.interp_x{k-1}, Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = TTFun.local_error(obj.data.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = TTFun.local_block(obj.base.oneds{k}, obj.data.res_x{k-1}, Jr_right, func);
                    obj.n_eval = obj.n_eval + nf + nr;
                    % different from left iteration, k + 1 is the previous index
                    % evaluate update function at x_k nodes
                    if k < d
                        [Fu,nu] = TTFun.local_block(obj.base.oneds{k}, obj.data.res_x{k-1}, Jx_right, func);
                        obj.n_eval = obj.n_eval + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.res_w{k}, obj.data.res_x{k}, obj.data.cores{k-1}] = TTFun.build_basis_amen(...
                        obj.base.oneds{k}, Jx_right, Jr_right, obj.data.res_w{k-1}, Rw_right, obj.data.cores{k-1}, F, Fu, Fr, ...
                        obj.data.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k-1) = size(obj.data.cores{k}, 1);
                end
            end
    end
    als_iter = als_iter + 1;
    %
    if als_iter == obj.opt.max_als || max(errs(ind)) < obj.opt.als_tol || obj.l2_err < obj.opt.als_tol
        % compute the last block
        %
        if obj.data.direction > 0 % k = d
            Jx_left = obj.data.interp_x{d-1};
            Jx_right = [];
            % evaluate interpolant function at x_k nodes
            [F, nf] = TTFun.local_block(obj.base.oneds{d}, Jx_left, [], func);
            % update relative error and increase counter
            %errs(d) = TTFun.local_error(obj.data.cores{d}, F);
            obj.n_eval = obj.n_eval + nf;
            %{
            % Schmidt operator and interpolation, k + 1 is the previous index
            % push couple matrix to the right
            [obj.data.cores{d}, obj.data.interp_x{d}, obj.data.cores{d-1}] = TTFun.build_basis_svd(obj.base.oneds{d},...
                Jx_right, obj.data.cores{d-1}, F, -obj.data.direction, obj.opt.int_method, ...
                obj.opt.local_tol, obj.opt.max_rank);
            %}
            nbleft  = size(F, 1);
            nnodes  = size(F, 2);
            nbright = size(F, 3);
            obj.data.cores{d} = reshape(F, nbleft, nnodes, nbright, []); % nbright = 1;
        else
            Jx_left = [];
            Jx_right = obj.data.interp_x{2};
            % evaluate interpolant function at x_k nodes
            [F, nf] = TTFun.local_block(obj.base.oneds{1}, Jx_left, Jx_right, func);
            % update relative error and increase counter
            %errs(1) = TTFun.local_error(obj.data.cores{1}, F);
            obj.n_eval = obj.n_eval + nf;
            %{
            % Schmidt operator and interpolation, k - 1 is the previous index
            % push couple matrix to the right
            [obj.data.cores{1}, obj.data.interp_x{1}, obj.data.cores{2}] = TTFun.build_basis_svd(obj.base.oneds{1},...
                Jx_left, obj.data.cores{2}, F, -obj.data.direction, obj.opt.int_method, ...
                obj.opt.local_tol, obj.opt.max_rank);
            %}
            nbleft  = size(F, 1);
            nnodes  = size(F, 2);
            nbright = size(F, 3);
            obj.data.cores{1} = reshape(F, nbleft, nnodes, nbright, []); % nbleft = 1;
        end
        %
        % evaluate the TTFun and give error estimates
        [obj.l2_err,obj.linf_err] = rel_error(obj);
        % Print the information of each TT iteration
        if isdebug(obj.var)
            fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10.5E   %10.5E  %10i\n', ...
                [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), obj.linf_err, obj.l2_err, obj.n_eval]);
        else
            fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10i\n', ...
                [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), obj.n_eval]);
        end
        %
        disp('  >> ALS completed, TT ranks')
        disp(['  >> ' num2str(rs)])
        break;
    else
        % evaluate the TTFun and give error estimates
        [obj.l2_err,obj.linf_err] = rel_error(obj);
        % Print the information of each TT iteration
        if isdebug(obj.var)
            fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10.5E   %10.5E  %10i\n', ...
                [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), obj.linf_err, obj.l2_err, obj.n_eval]);
        else
            fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10i\n', ...
                [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), obj.n_eval]);
        end
        % flip direction
        obj.data.direction = -obj.data.direction;
    end
end

end