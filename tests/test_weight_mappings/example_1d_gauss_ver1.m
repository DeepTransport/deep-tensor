% one dim gaussian

fun = @(x) 0.5*x.^2;
g_fun = @(x) x;

scales = 1:9;

for i = 1:numel(scales)
    map_A{i} = AlgebraicMapping(scales(i));
    map_L{i} = LogarithmicMapping(scales(i));
end
c1 = Chebyshev1st(40);
le = Legendre(40);

xx = linspace(-5,5,1E4);

figure
f = fun(xx);
g = g_fun(xx);
for i = 1:numel(scales)
    [z,dzdx] = domain2reference(map_A{i}, xx);
    [dlogzdx,d2logzdx2] = domain2reference_log_density(map_A{i}, xx);
    %
    logw_C = eval_log_measure(c1, z);
    logw_L = eval_log_measure(le, z);
    log_mu_AC = logw_C + dlogzdx;
    log_mu_AL = logw_L + dlogzdx;
    %
    g_logw_C = eval_log_measure_deri(c1, z);
    g_logw_L = eval_log_measure_deri(le, z);
    g_AC = g_logw_C.*dzdx + d2logzdx2;
    g_AL = g_logw_L.*dzdx + d2logzdx2;
    
    figure(1)
    subplot(3,3,i)
    plot(xx, exp(- f - log_mu_AC) )
    hold on
    plot(xx, exp(- f - log_mu_AL) )
    title(num2str(scales(i)))
    
    figure(2)
    subplot(3,3,i)
    plot(xx, exp(-f-dlogzdx).*(g+d2logzdx2).^2 )
    title([num2str(scales(i)), ' ', num2str(sum(exp(-f-dlogzdx).*(g+d2logzdx2).^2)/length(xx))])
end

for i = 1:numel(scales)
    [z,dzdx] = domain2reference(map_A{i}, xx);
    [dlogxdz,d2logxdz2] = reference2domain_log_density(map_A{i}, z);
    %
    g_logw_C = eval_log_measure_deri(c1, z);
    g_logw_L = eval_log_measure_deri(le, z);
    g_AC = -g./dzdx + d2logxdz2 - g_logw_C;
    g_AL = -g./dzdx + d2logzdx2 - g_logw_L;
    
    figure(3)
    subplot(3,3,i)
    plot(xx, exp(-f).*(g_AC).^2 )
    title([num2str(scales(i)), ' ', num2str(sum(exp(-f).*(g_AC).^2)/length(xx))])
end

%%
zz = linspace(-1,1,1E4);
for i = 1:numel(scales)
    [x,dxdz] = reference2domain(map_A{i}, zz);
    [dlogxdz,d2logxdz2] = reference2domain_log_density(map_A{i}, zz);
    %
    logw_C = eval_log_measure(c1, zz);
    logw_L = eval_log_measure(le, zz);
    log_mu_AC = logw_C - dlogxdz;
    log_mu_AL = logw_L - dlogxdz;
    %
    g_logw_C = eval_log_measure_deri(c1, zz);
    g_logw_L = eval_log_measure_deri(le, zz);
    g_AC = g_logw_C - d2logxdz2;
    g_AL = g_logw_L - d2logxdz2;
    %
    f = fun(x);
    g = g_fun(x);
    g = g.*dxdz;
    
    figure(4)
    subplot(3,3,i)
    plot(zz, exp(- f - log_mu_AC) )
    hold on
    plot(zz, exp(- f - log_mu_AL) )
    title(num2str(scales(i)))
    
    figure(5)
    subplot(3,3,i)
    plot(xx, exp(-f+dlogxdz).*(g+d2logxdz2).^2 )
    title([num2str(scales(i)), ' ', num2str(sum(exp(-f+dlogxdz).*(g+d2logxdz2).^2)/length(xx))])
end

%%
%{
xx = linspace(-5,5,1E4);
f = fun(xx);
for i = 1:numel(scales)
    [z,dzdx] = domain2reference(map_L{i}, xx);
    %
    logw_C = eval_log_measure(c1, z);
    logw_L = eval_log_measure(le, z);
    log_mu_LC = logw_C + log(dzdx);
    log_mu_LL = logw_L + log(dzdx);
    
    figure(3)
    subplot(3,3,i)
    plot(xx, exp(- f - log_mu_LC) )
    hold on
    plot(xx, exp(- f - log_mu_LL) )
    title(num2str(scales(i)))
end

zz = linspace(-1,1,1E4);
for i = 1:numel(scales)
    [x,dxdz] = reference2domain(map_L{i}, zz);
    %
    logw_C = eval_log_measure(c1, zz);
    logw_L = eval_log_measure(le, zz);
    log_mu_LC = logw_C - log(dxdz);
    log_mu_LL = logw_L - log(dxdz);
    %
    f = fun(x);
    
    figure(4)
    subplot(3,3,i)
    plot(zz, exp(- f - log_mu_LC) )
    hold on
    plot(zz, exp(- f - log_mu_LL) )
    title(num2str(scales(i)))
end
%}