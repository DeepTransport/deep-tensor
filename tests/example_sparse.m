
% setup the reference
tol = 1E-10;
dom = BoundedDomain([tol,1]);
d = 2;

base = ApproxBases(Legendre(30), dom, d);


debug_size = 1E4;
debug_x = sample_measure(base, debug_size);
deb = InputData([],debug_x);

%%%%

func1 = @(x) sqrt(1./sum(1E-5+x.^2,1));
func2 = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
func3 = @(x) [(1 + sum(x,1)).^(-d-1); exp( - sum(abs(x - 0.5), 1)); cos( sum(x.^2,1) )];

func = func1;

%%%%

opt1 = SparseOption();
opt1.tol = 1e-2;
opt1.init_total_degree = 20;
%opt.indexset = 'hyperbolic';
opt1.max_dim_basis = 1e3;
opt1.max_sample_size = 1e4;
opt1.enrich_degree = 1;
opt1.init_sample_size = 4;
opt1.enrich_sample_size = 1; 
opt1.adaptation_rule = 'reducedmargin';

opt2 = opt1;
opt2.fast = true;

opt3 = opt1;
opt3.weight_rule = 'none';

f = @(z) ApproxFun.feval_reference(base, func, z);

spfun{1} = SparseFun(f, base, opt1, 'var', deb);
spfun{2} = SparseFun(f, base, opt2, 'var', deb);
spfun{3} = SparseFun(f, base, opt3, 'var', deb);
%
spfun{4} = SparseFun(f, spfun{1}, 'var', deb);

tic; exact = func(debug_x); toc
for i = 1:4
    tic; app{i} = eval(spfun{i}, debug_x); toc
end

for i = 1:4
    figure
    plot(exact(:) - app{i}(:), '.')
end

x = debug_x(:,10);
for i = 1:3
    [gx,fx] = grad(spfun{i}, x);
    f = eval(spfun{i}, x);
    tol = 1E-5;
    fp = zeros(size(gx));
    fm = zeros(size(gx));
    for ii = 1:length(x)
        xp = x;
        xp(ii) = xp(ii)+tol;
        xm = x;
        xm(ii) = xm(ii)-tol;
        fp(ii) = eval(spfun{i}, xp);
        fm(ii) = eval(spfun{i}, xm);
    end
    disp(norm((fp-fm)/(2*tol) - gx))
end


