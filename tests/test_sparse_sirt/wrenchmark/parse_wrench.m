% Parse Wrench parameters
function [params] = parse_wrench(varargin)
% Parse model parameters
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end

if (~isfield(params, 'd_theta'))
    params.d_theta = input('State reduced dimension d_theta = ? (default 7): ');
    if (isempty(params.d_theta))
        params.d_theta = 7;
    end
end
if (~isfield(params, 'd_y'))
    params.d_y = input('Data reduced dimension d_y = ? (default 26): ');
    if (isempty(params.d_y))
        params.d_y = 26;
    end
end
if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Standard deviation of the observation noise sigma_n = ? (default 0.1): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 0.1; % Noise std
    end
end

if (~isfield(params, 'Nsamples'))
    params.Nsamples = input('Number of posterior MCMC samples Nsamples = ? (default 1e4): ');
    if (isempty(params.Nsamples))
        params.Nsamples = 1e4;
    end
end

% Grid for DIRT Level 0
if (~isfield(params, 'n'))
    params.n = input('Polynomial degree n = ? (default 20): ');
    if (isempty(params.n))
        params.n = 20;
    end
end
% Tempering powers
if (~isfield(params, 'beta0'))
    params.beta0 = input('Initial tempering power beta0 = ? (default 1e-3): ');
    if (isempty(params.beta0))
        params.beta0 = 1e-3;
    end
end
if (~isfield(params, 'ess_tol'))
    params.ess_tol = input('Tolerance for ESS between levels ess_tol = ? (default 0.5): ');
    if (isempty(params.ess_tol))
        params.ess_tol = 0.5;
    end
end
end