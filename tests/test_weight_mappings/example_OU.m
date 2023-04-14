
% setup the OU process
d = 10;
a = 0.8;

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

scales = 1:8;

for i = 1:numel(scales)
    doms_A{i} = AlgebraicMapping(scales(i));
    doms_L{i} = LogarithmicMapping(scales(i));
    %
    bases{i,1} = ApproxBases(Chebyshev1st(40), doms_A{i}, d);
    bases{i,2} = ApproxBases(Chebyshev1st(40), doms_L{i}, d);
    %
    bases{i,3} = ApproxBases(Legendre(40), doms_A{i}, d);
end

for i = 1:numel(scales)
    for j = 1:3
        %
        [z,dzdx] = domain2reference(bases{i,j}, deb.samples);
        mlogw = eval_measure_potential_reference(bases{i,j}, z);
        %
        pot = deb.f - mlogw + sum(log(dzdx),1);
        %
        vp_b(i,j) = var(pot);
        mp_b(i,j) = mean(pot);
        ve_b(i,j) = var(exp(-0.5*pot));
        me_b(i,j) = mean(exp(-0.5*pot));
    end
end


for i = 1:numel(scales)
    for j = 1:3
        %
        [z, mlogf] = sample_measure_reference(bases{i,j},1E4);
        [x, dxdz] = reference2domain(bases{i,j},z);
        y = feval(func, x);
        %
        pot = y - sum(log(dxdz),1);
        %
        vp_f(i,j) = var(pot);
        mp_f(i,j) = mean(pot);
        ve_f(i,j) = var(exp(-0.5*pot));
        me_f(i,j) = mean(exp(-0.5*pot));
    end
end

%return

opt = FTTOption('tt_method', 'random', 'als_tol', 1E-4, 'local_tol', 1E-20, ...
    'init_rank', 20, 'max_rank', 20, 'kick_rank', 0, 'max_als', 2);

for i = 1:numel(scales)
    disp(i)
    for j = 1:3
        tic;
        irts{i,j} = TTSIRT(func, bases{i,j}, opt, 'debug', deb, 'sample_x', sample_x);
        toc
    end
end


