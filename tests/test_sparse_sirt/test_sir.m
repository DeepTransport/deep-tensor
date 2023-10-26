function test_sir
addpath('err_estimate');
d = 2; % number of compartments
nruns = 1; % number of parallel runs

% Adjacency matrix of the diffusion
W = spdiags(ones(d,1)*[1 -2 1], -1:1, d, d);
W(1,d) = 1;
W(d,1) = 1;
W = W*0.5;
if (d==2)
    W = [-1 1; 1 -1];
end
if (d==1)
    W = 0;
end

d = size(W,1);

% True parameters
beta = 0.1*ones(1,d);
gamma = 1*ones(1,d);

% Try different initial states in different compartments
x0 = [99; 0; 0] + [0;1;0]*(d:-1:1);
x0(1,:) = 100-x0(2,:);
x0 = x0(:);

% Observation times
tobs = linspace(0,5,7);
% True data
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tobs, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));

% % fine discretization we need for events
% tfine = linspace(0,5,25*6+1);
% ind_obs = 25+1:25:25*6+1;
% assert(norm(tfine(ind_obs)-tobs(2:end))==0)

% plot(t,x);
sigma_n = 1;
data = x(2:end, 2+(0:d-1)*3);
% data = data + sigma_n*randn(size(data));
data = data(:);

func = @(bg)sir_ll_ftt(bg, data, W, sigma_n, x0, tobs);

% betas = 10.^(-3:1:0);
% temp1 = Tempering1(betas);

temp1 = Tempering1('min_beta', 1E-3, 'ess_tol', 0.1, 'ess_tol_init', 0.8);

% temp1 = DataBatchAdapt(1:numel(data), 'ess_tol', 0.1, 'ess_tol_init', 0.8);
% temp1 = DataBatchFixed({3, [1:2], [3:4], [5:6]}, 'ess_tol', 0.1, 'ess_tol_init', 0.8);

p = 30;
dom = BoundedDomain([0,2]);
base = ApproxBases(Legendre(p), dom, 2*d);
diag = UniformReference(dom);
%diag = GaussReference(0, 1, dom);

dirt_opt = DIRTOption('method', 'Eratio', 'defensive', 1E-13);


% %% My test
% test_tt_dirt(d, @(bgm)sir_ll(bgm, data, W, sigma_n, x0, tobs), betas)

%% Poly-DIRT approx
opt = SparseOption();
opt.tol = 1e-1;
opt.init_total_degree = 5;
%opt.indexset = 'hyperbolic';
opt.max_dim_basis = 2e4;
opt.max_sample_size = 2e5;
opt.enrich_degree = 1;
opt.init_sample_size = 2;
opt.enrich_sample_size = 1; 
opt.fast = true; 
opt.adaptation_rule = 'margin';
opt.opt_tol = 10;
opt.bulk_parameter = 0.1;

%dirt_opt2 = DIRTOption('method', 'Eratio', 'defensive', 1E-13);
% dirt_opt2 = DIRTOption('method', 'Eratio', 'defensive', 1E-13, 'n_debugs', 0);
% init_samples = sample_measure(base, dirt_opt2.n_samples);

for irun=1:nruns
    tic;
    poly_dirt = SparseDIRT(func, base, temp1, diag, opt, dirt_opt); % , 'init_samples', init_samples);
    ttimes_poly(irun) = toc;
    
    z = random(diag, 2*d, 2^13);
    mlogf = zeros(1, size(z,2));
    for i=1:size(z,2)
        [z(:,i), mlogf(i)] = eval_irt(poly_dirt, z(:,i));
    end
    [mllkd, mlp] = func(z); % , 1:numel(data));
    
%     z2 = mcmc_prune(z', -(mllkd+mlp)', -mlogf');
%     tau_poly(irun) = statsiact(z2)
    ess_poly(irun) = 1/ess_ratio(-(mllkd+mlp)' + mlogf')
    [~, hell_poly(irun), ~] = f_divergence(-mlogf, -(mllkd+mlp)); % hellinger(-(mllkd+mlp)', -mlogf')
    hell_poly(irun) = sqrt(hell_poly(irun))
    
%     scatter3(z(1,:), z(2,:), exp(-mllkd-mlp+mlogf), 1, exp(-mllkd-mlp+mlogf)); view(2)
    
    evals_poly(irun) = poly_dirt.n_eval
    Nbetas_poly(irun) = length(poly_dirt.bridge.betas);
end


%% FTT-DIRT Approx
% init_samples = sample_measure(base, dirt_opt1.n_samples+dirt_opt1.n_debugs);


for irun=1:nruns
    r = ttrankest(base, evals_poly(irun)/Nbetas_poly(irun)/2)
    tt_opt = TTOption('tt_method', 'fix_rank', ...
        'als_tol', 1E-4, 'local_tol', 1E-10, 'kick_rank', 0, 'max_rank', r, 'max_als', 2, 'init_rank', r);

    tic;
    tt_dirt = TTDIRT(func, base, temp1, diag, tt_opt, dirt_opt); % , 'init_samples', init_samples);
    ttimes_tt(irun) = toc;

    [z,mlogf] = eval_irt(tt_dirt, random(diag, 2*d, 2^13));
    [mllkd, mlp] = func(z); % , 1:numel(data));
    
%     z2 = mcmc_prune(z', -(mllkd+mlp)', -mlogf');
%     tau_tt(irun) = statsiact(z2)
    ess_tt(irun) = 1/ess_ratio(-(mllkd+mlp)'  + mlogf')
    [~, hell_tt(irun), ~] = f_divergence(-mlogf, -(mllkd+mlp)); % hellinger(-(mllkd+mlp)', -mlogf')    
    hell_tt(irun) = sqrt(hell_tt(irun))
    
%     scatter3(z(1,:), z(2,:), exp(-mllkd-mlp+mlogf), 1, exp(-mllkd-mlp+mlogf)); view(2)
    
    evals_tt(irun) = tt_dirt.n_eval   
    Nbetas_tt(irun) = length(tt_dirt.bridge.betas);
end

fprintf(['\t\t%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g \n'], ...
        mean(Nbetas_tt), std(Nbetas_tt), ...
        mean(ttimes_tt), std(ttimes_tt), mean(evals_tt), std(evals_tt), ...
        mean(ess_tt), std(ess_tt), mean(hell_tt), std(hell_tt));
    
fprintf(['\t\t%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g \n'], ...
        mean(Nbetas_poly), std(Nbetas_poly), ...
        mean(ttimes_poly), std(ttimes_poly), mean(evals_poly), std(evals_poly), ...
        mean(ess_poly), std(ess_poly), mean(hell_poly), std(hell_poly));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end


function test_tt_dirt(d, func, betas)
xsf = repmat({linspace(0,2,30)'}, 2*d, 1);

tic;
IRT = tt_dirt_approx(xsf, @(x,b1,b2)func(x)*(b2-b1), ...
                     betas, 'nswp', 2, 'kickrank', 2, 'y0', 5, 'IRTdenom', true, ...
                     'reference', 'n3', 'boundary', true, 'interpolation', 'f', 'stoptol', 0.1, 'trunctol', 1e-2);
                                  % ^ n4 needed for antithetic sampling       
toc;

evals = IRT.evalcnt

q = randref(IRT.reference, 2^14, 2*d);
        
tic;
[z,lFapp,lFex] = tt_dirt_sample(IRT, q, func); 
toc;

scatter3(z(:,1), z(:,2), exp(lFex-lFapp), 1, exp(lFex-lFapp)); view(2); drawnow;
% Normalising constant
Zpost = mean(exp(lFex-lFapp))
% Error estimate of the posterior
z2 = mcmc_prune(z, lFex, lFapp);
% tau_post = statsiact(z2)
ess_post= essinv(lFex,lFapp)
hell_post = hellinger(lFex,lFapp)
end
