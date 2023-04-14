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
    [dlogxdz,d2logxdz2] = reference2domain_log_density(map_A{i}, z);
    %
    logw_C = eval_log_measure(c1, z);
    logw_L = eval_log_measure(le, z);
    %
    g_logw_C = eval_log_measure_deri(c1, z);
    g_logw_L = eval_log_measure_deri(le, z);
    g_AC = -g./dzdx + d2logxdz2 - g_logw_C;
    g_AL = -g./dzdx + d2logxdz2 - g_logw_L;
    
    figure(1)
    subplot(3,3,i)
    plot(xx, exp(- f + dlogxdz - logw_C) )
    hold on
    plot(xx, exp(- f + dlogxdz - logw_L) )
    title(num2str(scales(i)))
    
    figure(2)
    subplot(3,3,i)
    plot(xx, exp(-f).*(g_AC).^2 )
    title([num2str(scales(i)), ' ', num2str(sum(exp(-f).*(g_AC).^2)*(10/length(xx)))])
    
    figure(3)
    subplot(3,3,i)
    plot(z, exp(- f + dlogxdz - logw_C) )
    hold on
    plot(z, exp(- f + dlogxdz - logw_L) )
    title(num2str(scales(i)))
    set(gca, 'xlim', [-1,1])
    
    figure(4)
    subplot(3,3,i)
    plot(z, exp(-f).*(g_AC).^2 )
    set(gca, 'xlim', [-1,1])
end


zz = linspace(-1,1,1E4);
for i = 1:numel(scales)
    [x,dxdz] = reference2domain(map_A{i}, zz);
    [dlogxdz,d2logxdz2] = reference2domain_log_density(map_A{i}, zz);
    %
    logw_C = eval_log_measure(c1, zz);
    logw_L = eval_log_measure(le, zz);
    %
    g_logw_C = eval_log_measure_deri(c1, zz);
    g_logw_L = eval_log_measure_deri(le, zz);
    %
    f = fun(x);
    g = g_fun(x);
    g_AC = -g.*dxdz + d2logxdz2 - g_logw_C;
    g_AL = -g.*dxdz + d2logxdz2 - g_logw_L;
    
    figure(5)
    subplot(3,3,i)
    plot(zz, exp(- f + dlogxdz - logw_C) )
    hold on
    plot(zz, exp(- f + dlogxdz - logw_L) )
    title(num2str(scales(i)))
    
    figure(6)
    subplot(3,3,i)
    plot(zz, exp(-f+dlogxdz).*(g_AC).^2 )
    title([num2str(scales(i)), ' ', num2str(sum(exp(-f+dlogxdz).*(g_AC).^2)*(2/length(zz)))])
end
