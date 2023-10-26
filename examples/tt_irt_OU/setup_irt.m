
% setup the OU process
d = 20;
a = 0.8;

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

dom = BoundedDomain([-5,5]);
base = ApproxBases(Lagrange1(160), dom, d);

opts{1} = TTOption('tt_method', 'amen', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);
opts{2} = TTOption('tt_method', 'random', ...
    'als_tol', 1E-4, 'local_tol', 1E-10, 'max_rank', 19, 'max_als', 4);


irts = cell(1,3);
for j = 1:2
    tic;
    deb.count.i = 0;
    irts{j} = TTIRT(func, base, opts{j}, 'var', deb);
    toc
end
tmp = round(irts{1}.approx, 1E-2);
deb.count.i = 0;
irts{3} = TTIRT(func, tmp, 'var', deb);

for j = 1:3
    tic;
    irts{j} = marginalise(irts{j}, -1);
    toc
end