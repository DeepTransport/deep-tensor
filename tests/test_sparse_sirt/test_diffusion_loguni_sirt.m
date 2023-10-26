% Parse parameters of the inverse model
addpath('diffusion');
addpath('err_estimate');
params.npi = 17;
params.delta = 0.1;
params.sigma_n = 1e10;  % 1e-2
params.m0 = 7;  % 15
params.log2N = 13;
params.rpi = 5; % TT ranks in DIRT  % 20 for channel

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


%%
p = 20;
% dom = BoundedDomain([-sqrt(3),sqrt(3)]);
%dom = BoundedDomain([-5,5]);
% dom = LogarithmicMapping(4);
%base = ApproxBases(Legendre(p), dom, d);

dom = AlgebraicMapping(4);
base = ApproxBases(Chebyshev1st(p), dom, d);

r = 2;
beta = 1E-2;
func = @(theta)diffusion_lpfun_ftt(theta, phil, bound, W2, Mass_summed, Q_obs, beta*params.sigma_n, mask);

TTOpt = TTOption('tt_method', 'fix_rank', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', r, 'max_als', 2, 'init_rank', r, 'kick_rank', 0);

irt = TTSIRT(func, base, TTOpt);