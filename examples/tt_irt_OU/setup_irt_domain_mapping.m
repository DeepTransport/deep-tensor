
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
bases{1} = ApproxBases(Lagrange1(80), doms{1}, d);

doms{2} = AlgebraicMapping(4);
bases{2} = ApproxBases(Lagrange1(80), doms{2}, d);

doms{3} = LogarithmicMapping(4);
bases{3} = ApproxBases(Lagrange1(80), doms{3}, d);

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