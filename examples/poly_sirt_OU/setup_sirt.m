
global neval__

% setup the OU process
d = 3;
a = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = OU(d, a);
func = @(x) eval_potential(data, x);

%%%%

debug_size = 1E4;
debug_x = data.B\randn(d, debug_size);
sample_x = data.B\randn(d, 1E3);
deb = InputData(sample_x, debug_x);

%%%%

% setup the reference polynomial

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

% bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{1} = ApproxBases(repmat({le},1,d), repmat({BoundedDomain([-4,4])},1,d));
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

opt = SparseOption();
opt.tol = 1e-2;
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

for i = 1:length(bases)
    neval__ = 0;
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'var', deb);
    toc
    fprintf('Real nevals = %d\n', neval__);
end


for i = 1:length(bases)
    r = ttrankest(bases{i}, irts{i}.approx.n_eval)
    tt_opt = TTOption('tt_method', 'fix_rank', ...
        'als_tol', 1E-4, 'local_tol', 1E-10, 'kick_rank', 0, 'max_rank', r, 'max_als', 1, 'init_rank', r);
    neval__ = 0;
    tic;
    tt_irts{i} = TTSIRT(func, bases{i}, tt_opt, 'var', deb);
    toc
    fprintf('Real nevals = %d\n', neval__);
end


function r = ttrankest(bases, N)
% Solve the quadratic to estimate the TT rank matching N evaluations
% (n(1)+n(d))*r + sum(n(2:d-1))*r^2 = N
% ar^2 + br - N = 0          a = sum(n(2:d-1)), b = n(1)+n(d)
% D = b^2 + 4*a*N
% r = (-b + sqrt(D))/(2*a)

a = sum(cellfun(@(b)numel(b.nodes), bases.oneds(2:end-1)));
b = numel(bases.oneds{1}.nodes) + numel(bases.oneds{end}.nodes);
r = (-b + sqrt(b^2 + 4*a*N))/(2*a);
r = round(r);
end

