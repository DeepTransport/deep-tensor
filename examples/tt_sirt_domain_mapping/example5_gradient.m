
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);

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

%{
doms{1} = BoundedDomain([-5,5]);
bases{1} = ApproxBases(Lagrangep(5,8), doms{1}, d);

doms{2} = AlgebraicMapping(4);
bases{2} = ApproxBases(Chebyshev1st(40), doms{2}, d);

doms{3} = AlgebraicMapping(4);
bases{3} = ApproxBases(Lagrangep(5,8), doms{3}, d);

doms{4} = LogarithmicMapping(4);
bases{4} = ApproxBases(Chebyshev1st(40), doms{4}, d);

doms{5} = LogarithmicMapping(4);
bases{5} = ApproxBases(Legendre(40), doms{5}, d);

opt = FTTOption('tt_method', 'random', 'als_tol', 1E-4, 'local_tol', 1E-20, ...
    'init_rank', 20, 'max_rank', 20, 'kick_rank', 0, 'max_als', 2);

for i = 1:length(bases)
    tic;
    irts{i} = TTSIRT(func, bases{i}, opt, 'debug', deb, 'samples', sample_x);
    toc
end

for i = 1:length(bases)
    tic;
    irts{i} = marginalise(irts{i}, -1);
    toc
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = rand(d, 1E2);

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