function obj = build_lis(obj, func, bases, sirt_opt)

obj.irts = {};
obj.n_eval = 0;
obj.logz = 0;

if ~isa(obj.ref, 'GaussReference')
    warning('Need a GaussReference for dimension reduction')
end
ns = get_pre_sample_size(obj);
samples = random(obj.ref, obj.d, ns);
        
n_layers = num_layers(obj.bridge); % should be 0
%
while n_layers < obj.dirt_opt.max_layers
    %
    [mllkds, mlps, mlogfx, gz_func, gz_ref, gz_irt] = ...
        eval_with_lis(obj.bridge, obj, func, obj.lis_opt.ratio_method, samples);
    if n_layers == 0
        obj.bridge = set_init(obj.bridge, mllkds);
    end
    %
    obj.bridge = adapt_density(obj.bridge, obj.dirt_opt.method, mllkds, mlps, mlogfx);
    %
    [str,log_weights] = print_str(obj.bridge, mllkds, mlps, mlogfx);
    %
    [V,s,~] = build_lis(obj.bridge, obj.lis_opt.ratio_method, gz_func, gz_ref, gz_irt, log_weights);
    %
    disp([str, sprintf(', #fevals=%3.3e \n', obj.n_eval)]);
    %
    switch obj.lis_opt.method
        case {'reduction'}
            cs = cumsum(s.^2);
            err = 0.5*sqrt(cs(end)-cs);
            r = min(sum(err>obj.lis_opt.tol), obj.lis_opt.max_rank);
            %
            %r = min(sum(s>obj.lis_opt.tol), obj.lis_opt.max_rank);
            r = max(r, obj.lis_opt.min_rank);
            %
            obj.basis{n_layers+1} = V(:,1:r);
            obj.basis_r(n_layers+1) = r;
            obj.basis_s{n_layers+1} = s;
            disp(['LIS dimension: ' num2str(r)])            
        case {'unitary'}
            if size(V,2) < size(V,1)
                warning('cannot form a unitary basis, use truncation instead');
            end
            obj.basis{n_layers+1} = V;
            obj.basis_r(n_layers+1) = size(V,2);
            obj.basis_s{n_layers+1} = s;
    end
    %
    ind = datasample(1:numel(log_weights), numel(log_weights), 'weights', exp(log_weights),'Replace',true);
    %
    %mlogfx = 0.5*sum((obj.basis{n_layers+1}'*samples).^2, 1);
    density = get_ratio_fun(obj.bridge, obj, obj.basis{n_layers+1}'*samples, mllkds, mlps, mlogfx);
    %
    obj.irts{n_layers+1} = get_new_layer(obj, func, bases, sirt_opt, n_layers, ...
        obj.basis{n_layers+1}'*samples(:,ind), density(ind));
    %
    obj.logz = obj.logz + log(obj.irts{n_layers+1}.z);
    obj.n_eval = obj.n_eval + obj.irts{n_layers+1}.approx.n_eval;
    
    % stop
    if islast(obj.bridge)
        break;
    end
    n_layers = num_layers(obj.bridge);
end

% evaluate the map and the target function
[x,mlogfx] = eval_irt(obj, samples);
if nargin(func) == 1
    [mllkds, mlps] = func(x);
elseif nargin(func) == 2
    [mllkds, mlps] = func(x, []);
end

% squared Hellinger error between logfx and (-mllkds-mlps)
[~,dh2,~] = f_divergence(-mlogfx, -mllkds-mlps);

%
fprintf('\n\niteration=%2d, Hell error=%3.3e, cum#fevals=%3.3e\n', ...
    n_layers, sqrt(dh2), obj.n_eval);
%
end