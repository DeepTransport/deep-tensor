function obj = build(obj, func, bases, sirt_opt)    

obj.irts = {};
obj.n_evals = 0;
obj.logz = 0;
%
ns = get_pre_sample_size(obj);

n_layers = num_layers(obj.bridge); % should be 0
%    
while n_layers < obj.dirt_opt.max_layers
    %
    if n_layers == 0
        % evaluate the map and the target function
        if isempty(obj.init_samples)
            [samples, mlogfx] = sample_measure(bases{1}, ns);
            [mllkds, mlps] = eval_wo_lis(obj.bridge, func, samples);
        else
            samples = obj.init_samples;
            [mllkds, mlps] = eval_wo_lis(obj.bridge, func, samples);
            mlogfx = mlps; % to counter balance the prior
            obj.bridge = set_init(obj.bridge, mllkds);
        end
    else
        % reuse samples
        [x,mlogfx] = eval_irt(obj, samples);
        [mllkds, mlps] = eval_wo_lis(obj.bridge, func, x);
    end
    %
    obj.bridge = adapt_density(obj.bridge, obj.dirt_opt.method, mllkds, mlps, mlogfx);
    %
    density = get_ratio_fun(obj.bridge, obj, samples, mllkds, mlps, mlogfx);
    %
    [str,log_weights] = print_str(obj.bridge, mllkds, mlps, mlogfx);
    %
    disp([str, sprintf(', #fevals=%3.3e \n', obj.n_evals)]);
    %
    ind = datasample(1:numel(log_weights), numel(log_weights), 'weights', exp(log_weights),'Replace',true);
    %
    % for the LIS
    obj.basis{n_layers+1} = [];
    obj.basis_r(n_layers+1) = 0;
    obj.basis_s{n_layers+1} = [];
    %
    obj.irts{n_layers+1} = get_new_layer(obj, func, bases, sirt_opt, n_layers, samples(:,ind), density(ind));
    %
    obj.logz = obj.logz + log(obj.irts{n_layers+1}.z);
    obj.n_evals = obj.n_evals + obj.irts{n_layers+1}.approx.n_evals;
    %
    % We need reference samples already here, in case we quit after 1 layer
    samples = random(obj.ref, obj.d, ns);
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
    num_layers(obj.bridge), sqrt(dh2), obj.n_evals);
%


end