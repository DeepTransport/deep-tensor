function scores = estimate_scaling(fun, bases, p, n)

xx = randn(2,n);
mlog_sampling = 0.5*sum(xx.^2,1);
%
[fxx, gxx] = fun(xx);
fxx = fxx - mlog_sampling;

scores = [];
for i = 1:numel(bases)
    [z,dzdx] = domain2reference(bases{i}, xx);
    %[dlogxdz,d2logxdz2] = reference2domain_log_density(bases{i}, z);
    %mlogw = eval_measure_potential_reference(bases{i}, z);
    %
    [~,d2logxdz2] = reference2domain_log_density(bases{i}, z);
    %
    g_mlogw = eval_measure_potential_reference_grad(bases{i}, z);
    
    g = -gxx./dzdx + d2logxdz2 + g_mlogw;
    scores = [scores, log(mean(exp(-fxx).*abs(g).^p, 2))/p];
end