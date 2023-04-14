
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

%{
% setup the reference polynomial

p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);

bases{1} = ApproxBases(le, BoundedDomain([-4,4]), d);
bases{2} = ApproxBases(c1, LogarithmicMapping(4), d);

%

opt = SparseOption();
opt.tol = 1e-2;
opt.init_total_degree = 3;
opt.max_dim_basis = 2e4;
opt.max_sample_size = 2e5;
opt.enrich_degree = 1;
opt.init_sample_size = 3;
opt.enrich_sample_size = 1; 
opt.adaptation_rule = 'margin';


for i = 1:length(bases)
    tic;
    irts{i} = SparseSIRT(func, bases{i}, opt, 'debug', deb);
    toc
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = rand(d, 123);

for i = 1:length(bases)
    [r,~,g] = eval_irt(irts{i}, z);
    tol = 1e-8;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        rp = r;
        rp(j,:) = rp(j,:) + tol;
        fp(j,:) = eval_potential(irts{i}, rp);
        
        rm = r;
        rm(j,:) = rm(j,:) - tol;
        fm(j,:) = eval_potential(irts{i}, rm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end

for i = 1:length(bases)
    [r,f,g] = eval_irt_reference(irts{i}, z);
    tol = 1e-8;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        rp = r;
        rp(j,:) = rp(j,:) + tol;
        fp(j,:) = eval_potential_reference(irts{i}, rp);
        
        rm = r;
        rm(j,:) = rm(j,:) - tol;
        fm(j,:) = eval_potential_reference(irts{i}, rm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end

i = 2;
r = eval_irt_reference(irts{i}, z);
g = eval_measure_potential_reference_grad(irts{i}.approx, r);
tol = 1e-8;
fp = zeros(size(g));
fm = zeros(size(g));
for j = 1:d
    rp = r;
    rp(j,:) = rp(j,:) + tol;
    fp(j,:) = eval_measure_potential_reference(irts{i}.approx, rp);
    
    rm = r;
    rm(j,:) = rm(j,:) - tol;
    fm(j,:) = eval_measure_potential_reference(irts{i}.approx, rm);
end
gd = ( fp - fm ) / (tol*2);
norm(g-gd,'fro')


%{
for i = 1:length(bases)
    [r,f,g] = eval_irt_reference(irts{i}, z);
    J = eval_rt_jac_reference(irts{i}, r, z);
    for j = 1:size(z,2)
        ind = (j-1)*d + (1:d);
        % J(:,ind) % dzdx, dfdx/dzdx
        g(:,j) = J(:,ind)'\g(:,j);
    end
    
    tol = 1e-4;
    fp = zeros(size(g));
    fm = zeros(size(g));
    for j = 1:d
        zp = z;
        zp(j,:) = zp(j,:) + tol;
        [~,fp(j,:)] = eval_irt_reference(irts{i}, zp);
        
        zm = z;
        zm(j,:) = zm(j,:) - tol;
        [~,fm(j,:)] = eval_irt_reference(irts{i}, zm);
    end
    gd = ( fp - fm ) / (tol*2);
    norm(g-gd,'fro')
end
%}