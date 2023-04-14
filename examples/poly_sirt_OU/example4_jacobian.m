
% setup the OU process
d = 3;
a = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = setup_ou_process(d, a);
func = @(x) eval_ou_process(data, x);

%%%%

debug_size = 1E4;
tmp = data.B\randn(d, debug_size);
tmp_f = feval(func, tmp);
deb = Debugger(tmp, tmp_f);
sample_x = data.B\randn(d, 1E3);

%%%%

% setup the reference polynomial

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

%
%{
opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 3;
opt.max_dim_basis = 2e3;
opt.max_sample_size = 2e4;
opt.enrich_degree = 1;
opt.init_sample_size = 3;
opt.enrich_sample_size = 1;
opt.adaptation_rule = 'reducedmargin';


for i = 1:length(bases)
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'debug', deb);
    toc
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = rand(d, 123);
for i = 1:length(bases)
    irts{i} = set_defensive(irts{i}, 5E-6);
    debug_jac(irts{i}, z, [3 1 2]);
    debug_jac(irts{i}, z, [1 3 2]);
end

