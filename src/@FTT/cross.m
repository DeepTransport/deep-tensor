function obj = cross(obj, func, sample_x, deb)
% Cross iterations for building the TT, only call this function if you know
% what you are doing

d = ndims(obj);
if obj.opt.kick_rank == 0
    obj.opt.tt_method = 'fix_rank';
end
% initilise TT cores
if isempty(obj.cores)
    % copy properties
    obj.direction = 1;
    % nested interpolation points
    obj.cores = cell(d,1);
    %
    % use the user provided sample set
    if size(sample_x, 2) < obj.opt.init_rank
        disp('Not enough number of samples to initialise ftt')
        sample_x = sample_measure_reference(obj, obj.opt.init_rank);
    end
    for k = d:-1:1
        obj.interp_x{k} = sample_x(k:d, 1:obj.opt.init_rank);
    end
    sample_x(:,1:obj.opt.init_rank) = [];
    %
    % determine m
    %y = feval_reference(obj,func,obj.interp_x{k}(:,1));
    y = func(obj.interp_x{k}(:,1));
    m = size(y,1);
    % interpolation weights
    obj.cores{1} = rand(1, cardinal(obj.oneds{1}), obj.opt.init_rank, m);
    for k = 2:d-1
        obj.cores{k} = rand(obj.opt.init_rank, cardinal(obj.oneds{k}), obj.opt.init_rank);
    end
    obj.cores{d} = rand(obj.opt.init_rank, cardinal(obj.oneds{d}), 1);
else
    obj.direction = -obj.direction;
end
% initialise the residual coordinates for AMEN
if use_amen(obj) && isempty(obj.res_x)
    if obj.direction > 0 % direction has already been flipped
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            sample_x = sample_measure_reference(obj, obj.opt.kick_rank);
        end
        % nested interpolation points for res
        for k = d:-1:1
            obj.res_x{k} = sample_x(k:d, 1:obj.opt.kick_rank);
        end
    else
        if size(sample_x, 2) < obj.opt.kick_rank
            disp('Not enough number of user provided samples to enrich ftt')
            sample_x = sample_measure_reference(obj, obj.opt.kick_rank);
        end
        % nested interpolation points for res, need to double check
        for k = d:-1:1
            obj.res_x{k} = sample_x(1:k, 1:obj.opt.kick_rank);
        end
    end
end
% reinitialise the residual blocks for AMEN
if use_amen(obj) && isempty(obj.res_w)
    if obj.direction < 0 % direction has already been flipped
        for k = 1:(d-1)
            obj.res_w{k} = ones(obj.opt.kick_rank, size(obj.cores{k},3));
        end
        obj.res_w{d} = ones(size(obj.cores{d},1), obj.opt.kick_rank);
    else
        obj.res_w{1} = ones(obj.opt.kick_rank, size(obj.cores{1},3));
        for k = 2:d
            obj.res_w{k} = ones(size(obj.cores{k},1), obj.opt.kick_rank);
        end
    end
end
% 
if isempty(deb)
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  cum#fevals\n');
else
    fprintf('\n\n  >> ALS  max_local_E  mean_local_E  max_r  max_debug_E  mean_debug_E  cum#fevals\n');
end

% start
l2_err = 0;
f_evals = 0;
als_iter = 0;
rs = ones(1,d);
while true % run ALS
    if obj.direction > 0
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
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.interp_x{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, Jx_left, obj.interp_x{k+1}, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf;
                    % Schmidt operator and interpolation, k - 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_left, obj.cores{k+1}, F, obj.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                    else
                        Jx_right = obj.interp_x{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, obj.interp_x{k-1}, Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf;
                    % Schmidt operator and interpolation, k + 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_right, obj.cores{k-1}, F, obj.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
        case {'random'}
            % at the head, update the random enrichment set
            if size(sample_x, 2) < obj.opt.kick_rank
                sample_x = [];
                enrich = sample_measure_reference(obj, obj.opt.kick_rank);
            else
                enrich = sample_x(:,1:obj.opt.kick_rank);
                sample_x(:,1:obj.opt.kick_rank) = [];
            end
            % start
            errs = zeros(1,d);
            for k = ind
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                    else
                        Jx_left = obj.interp_x{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, Jx_left, obj.interp_x{k+1}, func);
                    [Fe,ne] = local_block(obj, k, Jx_left, enrich(k+1:d,:),   func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf + ne;
                    % Schmidt operator and interpolation, k - 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_left, obj.cores{k+1}, cat(3,F,Fe), obj.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                    else
                        Jx_right = obj.interp_x{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, obj.interp_x{k-1}, Jx_right, func);
                    % different from left iteration, k + 1 is the previous index
                    [Fe,ne] = local_block(obj, k, enrich(1:k-1,:),   Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    f_evals = f_evals + nf + ne;
                    % Schmidt operator and interpolation, k + 1 is the previous index
                    % push couple matrix to the right
                    [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
                        Jx_right, obj.cores{k-1}, cat(1,F,Fe), obj.direction, obj.opt.int_method, ...
                        obj.opt.local_tol, obj.opt.max_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
        case {'amen'}
            % start
            errs = zeros(1,d);
            for k = ind
                if obj.direction > 0
                    if k == 1
                        Jx_left = [];
                        Jr_left = [];
                        Rw_left = 1;
                    else
                        Jx_left = obj.interp_x{k-1};
                        Jr_left = obj.res_x{k-1};
                        Rw_left = obj.res_w{k-1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, Jx_left, obj.interp_x{k+1}, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = local_block(obj, k, Jr_left, obj.res_x{k+1}, func);
                    f_evals = f_evals + nf + nr;
                    % evaluate update function at x_k nodes
                    if k > 1
                        [Fu,nu] = local_block(obj, k, Jx_left, obj.res_x{k+1}, func);
                        f_evals = f_evals + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k+1}] = FTT.build_basis_amen(...
                        obj.oneds{k}, Jx_left, Jr_left, Rw_left, obj.res_w{k+1}, obj.cores{k+1}, F, Fu, Fr, ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k) = size(obj.cores{k}, 3);
                else
                    if k == d
                        Jx_right = [];
                        Jr_right = [];
                        Rw_right = 1;
                    else
                        Jx_right = obj.interp_x{k+1};
                        Jr_right = obj.res_x{k+1};
                        Rw_right = obj.res_w{k+1};
                    end
                    % evaluate interpolant function at x_k nodes
                    [F, nf] = local_block(obj, k, obj.interp_x{k-1}, Jx_right, func);
                    % update relative error and increase counter
                    errs(k) = FTT.local_error(obj.cores{k}, F);
                    % evaluate res function at x_k nodes
                    [Fr,nr] = local_block(obj, k, obj.res_x{k-1}, Jr_right, func);
                    f_evals = f_evals + nf + nr;
                    % different from left iteration, k + 1 is the previous index
                    % evaluate update function at x_k nodes
                    if k < d
                        [Fu,nu] = local_block(obj, k, obj.res_x{k-1}, Jx_right, func);
                        f_evals = f_evals + nu;
                    else
                        Fu = Fr;
                    end
                    [obj.cores{k}, obj.interp_x{k}, obj.res_w{k}, obj.res_x{k}, obj.cores{k-1}] = FTT.build_basis_amen(...
                        obj.oneds{k}, Jx_right, Jr_right, obj.res_w{k-1}, Rw_right, obj.cores{k-1}, F, Fu, Fr, ...
                        obj.direction, obj.opt.int_method, obj.opt.local_tol, obj.opt.max_rank, obj.opt.kick_rank);
                    rs(k-1) = size(obj.cores{k}, 1);
                end
            end
    end
    als_iter = als_iter + 1;
    % evaluate the ftt and give error estimates
    debug_errs = [];
    if ~isempty(deb)
        approx = eval_reference(obj, deb.samples);
        l2_err = mean((deb.f(:) - approx(:)).^2).^0.5/mean(deb.f(:).^2).^0.5;
        debug_errs(1) = max(abs(deb.f(:) - approx(:)))  / max(abs(deb.f(:)));
        debug_errs(2) = l2_err;
    end
    % Print the information of each TT iteration
    if isempty(debug_errs)
        fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10i\n', ...
            [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), f_evals]);
    else
        fprintf('  >> %3i  %10.5E   %10.5E  %5i  %10.5E   %10.5E  %10i\n', ...
            [als_iter, max(errs(ind)), mean(errs(ind)), max(rs), debug_errs(1), debug_errs(2), f_evals]);
    end
    %
    if als_iter == obj.opt.max_als || max(errs(ind)) < obj.opt.als_tol || (~isempty(deb) && l2_err < obj.opt.als_tol)
        disp('  >> ALS completed, TT ranks')
        disp(['  >> ' num2str(rs)])
        break;
    else
        % flip direction
        obj.direction = -obj.direction;
    end
end
obj.l2_err = l2_err;
obj.n_evals = f_evals;
end