% Parse parameters of the inverse model
addpath('diffusion');
addpath('err_estimate');
params.npi = 17;
params.delta = 0.1;
params.sigma_n = 1e10;  % 1e-2
params.m0 = 7;  % 15
params.log2N = 13;
params.rpi = 5; % TT ranks in DIRT  % 20 for channel

beta = 10.^[-3:1/2:0];  % -6 for channel

params.sigma = 1;

% Forward model
params.runs = 1;
params.corr_length = 1; % 0.5;
params.nu = 2;
params.meshlevel = 2;
params.y0 = 0.5;          % default 1.5
                          % special values:
                          %     666 - high-contrast semi-annulus inside
                          %     777 - same + connected to x1=0 and x1=1

% Build the discretization and KLE
tol_kle = 1e-2;
[~,pm,bound,W1g,W1m,spind,Pua,phi,lambda,Mass, W2] = build_kle_eig(params.meshlevel, 'DN', params.nu, params.corr_length, tol_kle, params.m0);

Mass_summed = cellfun(@(M)sum(M,1), Mass, 'uni', 0);
Mass_summed = reshape(Mass_summed, [], 1);
Mass_summed = cell2mat(Mass_summed);

% weighted KLE components
d = numel(lambda);
phil = full(phi*spdiags(sqrt(lambda), 0, d, d));


if (exist(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0), 'file')>0)
    % Load the same observation for all experiments
    fprintf('Found Q_obs file for nu=%g, ml=%d, sn=%g, m0=%d, y0=%g\nRemove it to regenerate the observations\n', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0);
    load(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0));
else
    % Simulate some data
    fprintf('Generating new observation data\n');
    if (params.y0==666)
        C = ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2>0.4^2) ...
            & ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2<=0.5^2) ...
            & (pm(1,:)<=0.7) & (pm(2,:)>0.3);  % Semi-annulus inside the domain
        C = 0.01 + 100*double(C);  % 25
        C = log(C);
    elseif (params.y0==777)
        C = ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2>0.4^2) ...
            & ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2<=0.5^2) ...
            & (pm(1,:)<=0.7) & (pm(2,:)>0.3);
        C = C | ( (pm(1,:)<=0.35) & pm(2,:)>0.45 & pm(2,:) < 0.55 ) ...
            | ( (pm(1,:)>0.7) & pm(2,:)>0.72 & pm(2,:) < 0.78 );   % Channels to boundaries
        C = 0.01 + 100*double(C);  % 25
        C = log(C);
    else
        % y0 is the value of all random variables
        C = (sum(phil,2) * params.y0)';
    end
    [Q_obs, C, U] = diffusion_QoI(C, [], bound, W2, Mass_summed);
    Q_obs = Q_obs + randn(size(Q_obs)) * sqrt(params.sigma_n);
    save(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0), 'Q_obs');
end

% temp1 = Tempering1(1);
temp1 = Tempering1('min_beta', min(beta), 'ess_tol', 0.5, 'ess_tol_init', 0.5);  mask = [];
% temp1 = DataBatchAdapt(1:numel(Q_obs), 'ess_tol', 0.5, 'ess_tol_init', 0.5);

p = 20;
% dom = BoundedDomain([-sqrt(3),sqrt(3)]);
% dom = BoundedDomain([-5,5]);
dom = LogarithmicMapping(2);
%dom = AlgebraicMapping(4);
%base = ApproxBases(Legendre(p), dom, d);
base = ApproxBases(Chebyshev2nd(p), dom, d);
diag = UniformReference(dom);

dirt_opt = DIRTOption('method', 'Aratio', 'defensive', 1E-13);

opt = SparseOption();
opt.tol = 2e-2;
opt.init_total_degree = 5;
%opt.indexset = 'hyperbolic';
opt.max_dim_basis = 2e4;
opt.max_sample_size = 2e5;
opt.enrich_degree = 1;
opt.init_sample_size = 2;
opt.enrich_sample_size = 1; 
opt.fast = true; 
opt.adaptation_rule = 'reducedmargin';
opt.opt_tol = 10;
opt.bulk_parameter = 0.1;

% for irun=1:params.runs 
%     %% Poly DIRT       
%     tic;
%     poly_dirt = SparseDIRT(@(theta)diffusion_lpfun_ftt(theta, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n, mask), ...
%                            base, temp1, diag, opt, dirt_opt); % , 'init_samples', init_samples2);
%     ttimes_poly(irun) = toc;
%     
%     z = random(diag, d, 2^params.log2N);
%     mlogf = zeros(1, size(z,2));
%     for i=1:size(z,2)
%         [z(:,i), mlogf(i)] = eval_irt(poly_dirt, z(:,i));
%     end
%     [mllkd, mlp] = diffusion_lpfun_ftt(z, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n);
%     
%     z2 = mcmc_prune(z', -(mllkd+mlp)', -mlogf');
%     tau_poly(irun) = statsiact(z2)
%     ess_poly(irun) = essinv(-(mllkd+mlp)', -mlogf')
%     hell_poly(irun) = hellinger(-(mllkd+mlp)', -mlogf')        
%     evalcnt_poly(irun) = poly_dirt.n_eval
%     Nbetas_poly(irun) = length(poly_dirt.bridge.betas);
% end % irun


for irun=1:params.runs             
    %% FTT DIRT
%     r = ttrankest(base, evalcnt_poly(irun)/Nbetas_poly(irun)/2)
    r = 1;
    tt_opt = TTOption('tt_method', 'fix_rank', ...
        'als_tol', 1E-4, 'local_tol', 1E-10, 'kick_rank', 0, 'max_rank', r, 'max_als', 2, 'init_rank', r);

    tic;
    tt_dirt = TTDIRT(@(theta)diffusion_lpfun_ftt(theta, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n, mask), ...
                     base, temp1, diag, tt_opt, dirt_opt); % , 'init_samples', init_samples1);
    ttimes_tt(irun) = toc;
    
    z = random(diag, d, 2^params.log2N);
    mlogf = zeros(1, size(z,2));
    for i=1:size(z,2)
        [z(:,i), mlogf(i)] = eval_irt(tt_dirt, z(:,i));
    end
%     [z, mlogf] = eval_irt(tt_dirt, z);
    [mllkd, mlp] = diffusion_lpfun_ftt(z, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n);
    
    % normalising constant
    Z_post(irun) = mean(exp(-(mllkd+mlp)' + mlogf'));
    
    z2 = mcmc_prune(z', -(mllkd+mlp)', -mlogf');
    tau_tt(irun) = statsiact(z2)
    ess_tt(irun) = essinv(-(mllkd+mlp)', -mlogf')
    hell_tt(irun) = hellinger(-(mllkd+mlp)', -mlogf')
    evalcnt_tt(irun) = tt_dirt.n_eval
    Nbetas_tt(irun) = length(tt_dirt.bridge.betas);
end    

fprintf(['%3.3e+%3.3e\t%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g\t%g+%g \n'], ...
        mean(Z_post), std(Z_post), mean(Nbetas_tt), std(Nbetas_tt), ...
        mean(ttimes_tt), std(ttimes_tt), mean(evalcnt_tt), std(evalcnt_tt), ...
        mean(tau_tt), std(tau_tt), mean(ess_tt), std(ess_tt), mean(hell_tt), std(hell_tt));
    
% fprintf('%3.3e\t', tt_dirt.bridge.betas); fprintf('\n');

fprintf(['\t\t%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g\t%g+%g \n'], ...
        mean(Nbetas_poly), std(Nbetas_poly), ...
        mean(ttimes_poly), std(ttimes_poly), mean(evalcnt_poly), std(evalcnt_poly), ...
        mean(tau_poly), std(tau_poly), mean(ess_poly), std(ess_poly), mean(hell_poly), std(hell_poly));
    
% fprintf('%3.3e\t', poly_dirt.bridge.betas); fprintf('\n');    
