
% setup the OU process
d = 20;
a = 0.5;

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

doms{1} = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Lagrangep(5,8), doms{1}, d);

doms{2} = AlgebraicMapping(4);
bases{2} = ApproxBases(Chebyshev1st(40), doms{2}, d);

doms{3} = AlgebraicMapping(4);
bases{3} = ApproxBases(Lagrangep(5,8), doms{3}, d);

doms{4} = LogarithmicMapping(4);
bases{4} = ApproxBases(Chebyshev2nd(40), doms{4}, d);

doms{5} = LogarithmicMapping(4);
bases{5} = ApproxBases(Legendre(40), doms{5}, d);

opt = TTOption('tt_method', 'random', 'als_tol', 1E-4, 'local_tol', 1E-20, ...
    'init_rank', 20, 'max_rank', 20, 'kick_rank', 0, 'max_als', 2);

irts = cell(length(bases),1);
for i = 1:length(bases)
    tic;
    irts{i} = TTSIRT(func, bases{i}, opt, 'var', deb);
    toc
end

for i = 1:length(bases)
    tic;
    irts{i} = marginalise(irts{i}, -1);
    toc
end