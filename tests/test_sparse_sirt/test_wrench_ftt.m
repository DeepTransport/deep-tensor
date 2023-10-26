% A script to test Wrenchmark elasticity example of conditional DIRT,
% using FTT.
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on.
%   Model parameters:
%       d_theta: reduced dimension for Theta (state parameters)
%       d_y: reduced dimension for Y (data)
%       sigma_n: standard deviation of the observation noise (must be >0)
%   Test parameters:
%       Ndata: number of data samples to test
%       Nsamples: length of posterior MCMC chain to produce
%                 For RMSE estimation, Nsamples should be a multiple of 4
%   Approximation (DIRT) parameters:
%       n: number of grid points in each variable for Lagrange1 basis
%       R: TT rank (all TT ranks are set to R+1)
%       beta: a vector of tempering powers (must increase and end at 1)
%
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   DIRT: DIRT class (contains all TT decompositions)
%   wrench: an instance of the Wrench class defining the forward PDE
%   Q_obs: last data vector
%   t_ess: N/ESS for Importance Sampling
%   t_iact_2: IACT of unpreconditioned pCN chain
%   t_iact_3: IACT of negatively correlated preconditioned pCN
%   t_iact_4: IACT of uncorrelated preconditioned pCN chain
%   t_iact_5: IACT of unpreconditioned NUTS chain
%   t_iact_6: IACT of preconditioned NUTS
%   lE1, err_lE1: mean and std of posterior log E in Importance Sampling
%   lE2, err_lE2: mean and std of posterior log E in unpreconditioned pCN
%   lE3, err_lE3: mean and std of posterior log E in negative pCN
%   lE4, err_lE4: mean and std of posterior log E in uncorrelated pCN
%   lE5, err_lE5: mean and std of posterior log E in unpreconditioned NUTS
%   lE6, err_lE6: mean and std of posterior log E in preconditioned NUTS



function test_wrench_ftt(varargin)
addpath('wrenchmark');
addpath('err_estimate');

params = parse_wrench(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation of the wrench %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmax = 0.15; % mesh size: 0.5=coarse , 0.15=fine , 0.08 = ULTRA fine
precomputeEverything = 1; % you need it to be 1 in order to compute gradients, and 2 for parallel precomputation

wrench = Wrench(hmax,precomputeEverything);


% plot the geometry & the mesh
figure(4)
subplot(2,1,1)
wrench.plotGeometry();

subplot(2,1,2)
wrench.plotMesh();

%% The parameter is the "log of YoungModulus field"
% get a realization from it (Gaussian distrib)...
logE = wrench.logE(); 
% ... and plot it
figure(4);
clf;
wrench.plot(logE)

%% Compute and plot the solution
% logE = wrench.logE(ones(size(wrench.Sigma12,1),1));
logE = wrench.logE();
wrench.plotSol(logE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The observation and its gradient %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vertical displacement of the solution on the left of the wrench
% type = 1;
% Von Mises stress at some point
% type = 2;
% Displacement along the line
type = 3;

logE = wrench.logE();
[y,g,XObs,L] = wrench.evalObs(logE,type);
clf;
plot(XObs(1,:),y)

%% Norm on the data space (if type=3)

RX = wrench.M + wrench.K;
[~,~,~,L] = wrench.evalObs(wrench.logE(),3);
% L = wrench.B'*L; % already implemented in evalObs
RYinv = full(L'*(RX\L));
RY = inv(RYinv);

% Then, measure distances in Im(L) using the metric induced by 
%      <.,.> = (.^T)*RY*(.)

%% Matrix H for dimension reduction
K=500;
H = zeros(wrench.dimParam);
data = zeros(size(XObs,2),K);
for k=1:K
    logE = wrench.logE();
    [y,g] = wrench.evalObs(logE,type);
    data(:,k) = y;
    H = H + g'*RY*g;
end
H = H/K;
CovData = cov(data');


%% Compute the *informed* directions
[Us,Ds] = svd( wrench.Sigma12'*H*wrench.Sigma12);
Ds = diag(Ds);
Us = wrench.Sigma12*Us;

subplot(2,2,1)
semilogy(Ds)
title('eigenvalues')
for i=1:3
    subplot(2,2,i+1)
    wrench.plot(Us(:,i))
    title(['Mode ' num2str(i)])
end

%% Compute the *informative* directions
A = chol(RY);
[Ud,Dd] = svd(A*CovData*A');
Dd = diag(Dd);
% Ud = A\Ud;  % I need orthogonal Ud' to project full data into subspace
subplot(1,2,1)
semilogy(Dd)
legend('eigenvalues')
subplot(1,2,2)
plot(XObs(1,:),Ud(:,1:6))
legend(arrayfun(@(x) ['Mode ' num2str(x)],1:4,'UniformOutput',0))

%% Sample in LIS

d_y = params.d_y;
d_theta = params.d_theta;
Lprior = 3.5; % Prior domain will be [-Lprior,Lprior]^d

xi = randn(d_theta,1);
logE = wrench.logE(Us(:,1:d_theta)*xi);
clf;
wrench.plotSol(logE)

% %% Surrogate forward map
% 
% % Coefficient
% [x0,w] = lgwt(9,-Lprior,Lprior);
% coeffLogE = amen_cross_s([size(wrench.Sigma12,1); numel(x0)*ones(d_theta,1)], @(ind)(Us(ind(:,1),1:d_theta)*reshape(x0(ind(:,2:d_theta+1)), d_theta, 1)), 1e-8, 'vec', false, 'kickrank', 8);
% coeffLogE = round(coeffLogE, 1e-13);
% 
% coeffE = amen_cross_s({coeffLogE}, @(x)log(1+exp(x)), 1e-5, 'kickrank', 0.4);
% 
% % Test coeff approx
% theta = randref(sprintf('n%g', Lprior), 1e3, d_theta);
% err_e = zeros(size(theta,1),1);
% E2 = tt_sample_lagr(chunk(coeffE,2,d_theta+1), repmat({x0},d_theta,1), theta);
% E1 = squeeze(coeffE{1});
% for i=1:size(theta,1)
%     Eex = log(1+exp(wrench.logE(Us(:,1:d_theta)*theta(i,:)')));
%     err_e(i) = norm(E1*E2(i,:)' - Eex)/norm(Eex);
% end
% figure(1);
% plot(err_e);
% title('error in coefficient');
% drawnow;
% 
% %% Full solution surrogate
% [u,~,acrossevals] = als_cross_parametric(coeffE, @(C)wrench.evalFEMSystem(C), 1e-4, 'kickrank', 0, 'nswp', 1, 'random_init', 1e3);
% 
% % Test solution approx
% theta = randref(sprintf('n%g', Lprior), 1e3, d_theta);
% err_u = zeros(size(theta,1),1);
% U2 = tt_sample_lagr(chunk(u,2,d_theta+1), repmat({x0},d_theta,1), theta);
% U1 = squeeze(u{1});
% for i=1:size(theta,1)
%     [Uex,~,~] = wrench.evalFEMSystem(log(1+exp(wrench.logE(Us(:,1:d_theta)*theta(i,:)'))));
%     Uex = Uex{1};
%     err_u(i) = norm(U1*U2(i,:)' - Uex)/norm(Uex);
% end
% figure(2);
% plot(err_u);
% title('error in full solution');
% drawnow;
% 
% %% L is the observation operator
% [y,g,Xobs,L] = wrench.evalObs(wrench.logE(),type);
% 
% L = L / 1e3;  % The original units are scaled off, making setting SNR inconvenient
% 
% % Forward map surrogate
% G = chunk(u,2,d_theta+1);
% U1 = squeeze(u{1});
% U1 = Ud(:,1:d_y)'*(L'*U1);
% G = U1*G;
% 
% G = round(G, 1e-3);
% 
% % Test forward map approx
% theta = randref(sprintf('n%g', Lprior), 1e3, d_theta);
% err_obs = zeros(size(theta,1),1);
% Gz = tt_sample_lagr(G, repmat({x0},d_theta,1), theta);
% for i=1:size(theta,1)
%     y = wrench.evalObs(wrench.logE(Us(:,1:d_theta)*theta(i,:)'),type) / 1e3;
%     y = Ud(:,1:d_y)'*y;
%     err_obs(i) = norm(Gz(i,:)' - y)/norm(y);
% end
% figure(3);
% plot(err_obs);
% title('error in forward map');
% drawnow;


% %% Check likelihood approximation on LIS
% 
sigma_n = params.sigma_n; % std, |Q_obs| ~ 1 after rescaling (original max(|Q_obs|)~1e3)
% 
% xi_true = randn(d_theta,1);
% Q_obs = wrench.evalObs(wrench.logE(Us(:,1:d_theta)*xi_true), type) / 1e3;
% Q_obs = Q_obs + sigma_n*randn(size(Q_obs));
% Q_obs = Ud(:,1:d_y)'*Q_obs;
% 
% xi2 = xi_true + sigma_n*randn(d_theta,1);
% 
% mllkd = [];
% err_lkd = [];
% % miny = inf; % this was used to check that SNR ~ 1/sigma_n
% for k=K:-1:1
%     eta = randn(size(Us,2)-d_theta, 1);
%     y = wrench.evalObs(wrench.logE(Us(:,1:d_theta)*xi2 + Us(:,d_theta+1:end)*eta), type) / 1e3;
%     mllkd(k) = sum((A*(y - Ud(:,1:d_y)*Q_obs)).^2)/(2*sigma_n^2);
%     mllkdr = sum((A*(Ud(:,1:d_y)*Ud(:,1:d_y)'*y - Ud(:,1:d_y)*Q_obs)).^2)/(2*sigma_n^2);
%     err_lkd(k) = mllkd(k)./mllkdr - 1;    
% %     miny = min(miny, min(A*y))
% end
% 
% figure(4);
% plot(1:K, mllkd./mean(mllkd) - 1, 1:K, err_lkd);
% title('mllkd on LIS error')
% drawnow;

%% DIRT approx parameters

dom = BoundedDomain([-5, 5]);
% dom = AlgebraicMapping(4);
dom2 = BoundedDomain([0,1]);
base = ApproxBases(Chebyshev1st(params.n), dom, d_theta);
base2 = ApproxBases(Chebyshev1st(params.n), dom2, d_theta);
ref2 = UniformReference(dom2);

AU = (A*Ud(:,1:d_y))';  % chol(Mass) on Data LIS
     
% Data points
xi_obs = randn(d_theta,1);
if (d_theta==7)
    % Take 1 of the selected true parameters
    xi_obs = [-0.1752    0.3963   -0.0280    0.2508   -0.9472    0.5423    0.6350]';
end
% if (params.Ndata==-2)&&(d_theta==7)
%     xi_obs = [-1.4659   -1.1939   -0.9693    0.6767   -0.6462    0.1765   -0.1971]';    
%     params.Ndata = 1;
% end

% Data
if (exist(sprintf('Q_wrench_dtheta%d_dy%d_sigman%g.mat', d_theta, d_y, sigma_n), 'file')>0)
    % Load the same observation for all experiments
    fprintf('Found Q_wrench file for d_theta=%d, d_y=%d, sigma_n=%g\nRemove it to regenerate the observations\n', d_theta, d_y, sigma_n);
    load(sprintf('Q_wrench_dtheta%d_dy%d_sigman%g.mat', d_theta, d_y, sigma_n));
else
    Q_obs = wrench.evalObs(wrench.logE(Us(:,1:d_theta)*xi_obs), type) / 1e3;
    Q_obs = Q_obs + sigma_n*randn(size(Q_obs));
    save(sprintf('Q_wrench_dtheta%d_dy%d_sigman%g.mat', d_theta, d_y, sigma_n), 'Q_obs');
end
Q_obs = Ud(:,1:d_y)'*Q_obs; % Transform data into reduced space

temp1 = Tempering1('min_beta', params.beta0, 'ess_tol', params.ess_tol, 'ess_tol_init', 0.9);
    
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
opt.adaptation_rule = 'reducedmargin';
opt.opt_tol = 10;
opt.bulk_parameter = 0.1;

%%%%%%%%%%%%%%%%%%%%%%% Algebraic map ??????????????????????
%%%%%%%%%%%%%%%%%%%%%%% test_SparseDIRT2.m

dirt_opt = DIRTOption('method', 'Aratio', 'defensive', 1E-13);
%%%%%%%%%%%%%%%%% init_samples: just randn?
% init_samples1 = randn(d_theta, dirt_opt1.n_samples+dirt_opt1.n_debugs);
% sample_measure(base, dirt_opt1.n_samples+dirt_opt1.n_debugs);

% dirt_opt2 = DIRTOption('method', 'Aratio', 'defensive', 1E-13, 'n_debugs', 0);
% init_samples2 = randn(d_theta, dirt_opt2.n_samples);
% sample_measure(base, dirt_opt2.n_samples);
Nsamples = params.Nsamples;

% fun = @(xi)sep_likelihood(xi,G,AU,sigma_n,x0, Q_obs);
fun = @(xi)sep_full_lkd(xi,wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n);

%% Poly DIRT
for irun=1:1
    ttt = tic;
    poly_dirt = SparseDIRT(fun, {base base2}, temp1, ref2, opt, dirt_opt); % , 'init_samples', init_samples2);
    
    z = random(ref2, d_theta, Nsamples);
    mlogf = zeros(1, size(z,2));
    for i=1:size(z,2)
        [z(:,i), mlogf(i)] = eval_irt(poly_dirt, z(:,i));
    end
    [mllkds,mlps] = fun(z);
    ttimes_poly(irun) = toc(ttt);
    
    z2 = mcmc_prune(z', -(mllkds+mlps)', -mlogf');
    tau_poly(irun) = statsiact(z2)
    ess_poly(irun) = essinv(-(mllkds+mlps)', -mlogf')
    hell_poly(irun) = hellinger(-(mllkds+mlps)', -mlogf') 
    evalcnt_poly(irun) = poly_dirt.n_eval;
    Nbetas_poly(irun) = numel(poly_dirt.bridge.betas);
end


%% FTT DIRT
for irun=1:1
    r = ttrankest(base, evalcnt_poly(irun)/Nbetas_poly(irun)/2)
    tt_opt = TTOption('tt_method', 'fix_rank', ...
        'als_tol', 1E-4, 'local_tol', 1E-10, 'kick_rank', 0, 'max_rank', r, 'max_als', 2, 'init_rank', r);

    ttt = tic;
    tt_dirt = TTDIRT(fun, {base base2}, temp1, ref2, tt_opt, dirt_opt); % , 'init_samples', init_samples1);
        
    z = random(ref2, d_theta, Nsamples);
    mlogf = zeros(1, size(z,2));
    for i=1:size(z,2)
        [z(:,i), mlogf(i)] = eval_irt(tt_dirt, z(:,i));
    end
    % compute the minus log likelihood and minus log prior
    [mllkds,mlps] = fun(z);
    ttimes_tt(irun) = toc(ttt);
    
    z2 = mcmc_prune(z', -(mllkds+mlps)', -mlogf');
    tau_tt(irun) = statsiact(z2)
    ess_tt(irun) = essinv(-(mllkds+mlps)', -mlogf')
    hell_tt(irun) = hellinger(-(mllkds+mlps)', -mlogf')
    evalcnt_tt(irun) = tt_dirt.n_eval
    Nbetas_tt(irun) = numel(tt_dirt.bridge.betas);
end    

fprintf(['%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g\t%g+%g \n'], ...
        mean(Nbetas_tt), std(Nbetas_tt), ...
        mean(ttimes_tt), std(ttimes_tt), mean(evalcnt_tt), std(evalcnt_tt), ...
        mean(tau_tt), std(tau_tt), mean(ess_tt), std(ess_tt), mean(hell_tt), std(hell_tt));
    
fprintf(['\t\t%g+%g \n' ...
        '%g+%g\t%g+%g \n' ...
        '%g+%g\t%g+%g\t%g+%g \n'], ...
        mean(Nbetas_poly), std(Nbetas_poly), ...
        mean(ttimes_poly), std(ttimes_poly), mean(evalcnt_poly), std(evalcnt_poly), ...
        mean(tau_poly), std(tau_poly), mean(ess_poly), std(ess_poly), mean(hell_poly), std(hell_poly));
    

%% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end



%LIS space only                  
function [mlk,mlp,gmlk,gmlp] = sep_likelihood(xi,G,AU,sigma_n,x0, Q_obs)
d_y = size(AU, 1);
if (nargin>5) || (~isempty(Q_obs)) % This function is already conditional
    xi = [repmat(Q_obs, 1, size(xi,2)); xi];
end
d_theta = size(xi,1) - d_y;
xi = xi'; % from FTT
Gxi = tt_sample_lagr(G, repmat({x0},d_theta,1), xi(:, d_y+1:end));
Resid = (Gxi - xi(:,1:d_y))*AU; % size N x Ddata_full
mlk = sum(Resid.^2, 2)/(2*sigma_n^2);    
mlp = sum(xi(:, d_y+1:end).^2, 2)/2;
mlk = mlk';
mlp = mlp';
%%%%%%%%%%%%%%%
mlk = mlk + mlp;
mlp = mlp*0;
%%%%%%%%%%%%%%
if (nargout>2)
    % G is in Lagrange basis. Differentiate each separately
    GradGxi = zeros(size(xi,1), d_y, d_theta);
    for i=1:d_theta
        GradGxi(:,:,i) = tt_sample_lagr_diff(G, repmat({x0},d_theta,1), xi(:, d_y+1:end), i);
    end
    GradGxi = permute(GradGxi, [3,1,2]);
    GradGxi = reshape(GradGxi, [], d_y);
    GradGxi = GradGxi*AU;
    % tracemult over data dim
    GradGxi = reshape(GradGxi, d_theta, size(xi,1), []);
    GradGxi = permute(GradGxi, [3,1,2]); % size DData_full x d x N
    Resid = reshape(Resid.', 1, [], size(xi,1)); % size 1 x DData_full x N
    gmlk = tracemult(Resid, (1:size(xi,1))', GradGxi);  % size 1 x d x N
    gmlk = reshape(gmlk, d_theta, size(xi,1));
    gmlk = gmlk/(sigma_n^2);
    
    gmlp = xi(:, d_y+1:end)';
    %%%%%%%%%%%%%%%%%
    gmlk = gmlk + gmlp;
    gmlp = gmlp*0;
end
end


%LIS space only
function [mllkd,mlp] = log_target_pullback_pcn(irt,func,ry,z)

mlf = pullback(irt, func, [repmat(ry,1,size(z,2));z]);
%
mlp = 0.5*sum(z.^2,1);
mllkd = mlf - mlp;

end


function [mllkd,mlp]=pullback_full_pcn(xi,wrench,type,irt,qref_obs,Us,A,AU,d_theta,sigma_n)
% xi is a hybrid variable (u_r, theta_\perp), u_r in the reference of DIRT
u = xi(1:d_theta, :);
xperp = xi(d_theta+1:end, :);

% re-implement pullback here since it doesn't allow extra variables
z = [repmat(qref_obs,1,size(u,2));u];
[x,logft] = eval_irt(irt, z);
Q_obs = x(1:numel(qref_obs),1);  % recover Q_obs in the original space
x = x(end-d_theta+1:end, :);  % this is x_r in the LIS
% compute the minus log likelihood and minus log prior
[mllkds,mlps] = sep_full_lkd([x; xperp],wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n);
% compute the reference density at z
logfz = log_joint_pdf(irt.ref, z);
mlf = (mllkds+mlps) + logft - logfz;

%
mlp = 0.5*sum(xi.^2,1);
mllkd = mlf - mlp;

end


function [mlf,gmlf]=pullback_full_nuts(xi,wrench,type,irt,qref_obs,Us,A,AU,d_theta,sigma_n)
% xi is a hybrid variable (u_r, theta_\perp), u_r in the reference of DIRT
u = xi(1:d_theta, :);
xperp = xi(d_theta+1:end, :);

% re-implement pullback here since it doesn't allow extra variables
z = [repmat(qref_obs,1,size(u,2));u];
[x,logft,gz,Juz,Jux] = eval_irt(irt, z);
Q_obs = x(1:numel(qref_obs),1);  % recover Q_obs in the original space
x = x(end-d_theta+1:end, :);  % this is x_r in the LIS
% compute the minus log likelihood and minus log prior
[mllkds,mlps,gmllkds,gmlps] = sep_full_lkd([x; xperp],wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n);
% compute the reference density at z
[logfz,glfz] = log_joint_pdf(irt.ref, z);
% backtrack minus lpost and its gradient
mlf = (mllkds+mlps) + logft - logfz;
gmlf = DIRT.backtrack(Juz, Jux, [zeros(numel(qref_obs), size(u,2)); (gmllkds(1:d_theta,:)+gmlps(1:d_theta,:))]) + gz - glfz;
gmlf = [gmlf(end-d_theta+1:end, :); (gmllkds(d_theta+1:end,:)+gmlps(d_theta+1:end,:))];
end


function [mlf,gmlf] = lkd_for_nuts(xi,wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n)
if (nargin>1)
    [mlk,mlp,gmlk,gmlp]=sep_full_lkd(xi,wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n);
    gmlf = gmlk + gmlp;
else
    [mlk,mlp]=sep_full_lkd(xi,wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n);
end
mlf = mlk + mlp;
end

% Likelihood in the full space
function [mlk,mlp,gmlk,gmlp]=sep_full_lkd(xi,wrench,type,Q_obs,Us,A,AU,d_theta,sigma_n)
% xi is a hybrid variable (theta_r, theta_\perp), theta_r in LIS
xi = max(xi, -12);
xi = min(xi,  12);

u = xi(1:d_theta, :);

mlk = zeros(1,size(u,2));
mlp = zeros(1,size(u,2));
if (nargout>2)
    gmlk = zeros(size(Us,2), size(u,2));
    gmlp = xi; % [u; zeros(size(Us,2)-d, size(u,2))];
end
for k=1:size(u,2)
    thetafull = Us(:, 1:min(size(Us,2), size(xi,1))) * xi(:,k);
    if (nargout>2)
        [y,g] = wrench.evalObs(wrench.logE(thetafull), type);
        y = y / 1e3;
        g = g / 1e3;  % size of g is ddata x dtheta
        g = A*g;
    else
        y = wrench.evalObs(wrench.logE(thetafull), type) / 1e3;
    end
    Resid = A*y - AU'*Q_obs;
    mlk(k) = sum(Resid.^2)/(2*sigma_n^2);
    mlp(k) = sum(xi(:,k).^2, 1)/2;
    mlk(k) = mlk(k) + mlp(k); %%%%%%%%%%%%%%
    mlp(k) = 0; %%%%%%%%%%%
    if (nargout>2)
        gmlk(:,k) = (Us'*(g'*Resid))/(sigma_n^2);
        gmlk(:,k) = gmlk(:,k) + gmlp(:,k); %%%%%%%%%%%
        gmlp(:,k) = 0; %%%%%%%%%%%%%
    end
end
end



