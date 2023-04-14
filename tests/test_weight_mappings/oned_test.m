function oned_test(fun, g_fun, poly, maps, p)

% one dim gaussian

xx = linspace(-5,5,1E4);

figure
fxx = fun(xx);
gxx = g_fun(xx);
for i = 1:numel(maps)
    [z,dzdx] = domain2reference(maps{i}, xx);
    [dlogxdz,d2logxdz2] = reference2domain_log_density(maps{i}, z);
    %
    logw = eval_log_measure(poly, z);
    %
    g_logw = eval_log_measure_deri(poly, z);
    g = -gxx./dzdx + d2logxdz2 - g_logw;
    
    figure(1)
    subplot(3,3,i)
    plot(xx, exp(- fxx + dlogxdz - logw) )
    title(num2str(i))
    
    figure(2)
    subplot(3,3,i)
    plot(z, exp(- fxx + dlogxdz - logw) )
    title(num2str(i))
    set(gca, 'xlim', [-1,1])
    
    figure(3)
    subplot(3,3,i)
    plot(xx, exp(-fxx).*g.^2 )
    title([num2str(i), ' ', num2str(sum(exp(-fxx).*g.^2)*(10/length(xx)))])
    
    %pre_score = exp(-fxx+logw).*(-gxx./dzdx + d2logxdz2 - g_logw).^2;
    pre_score = exp(-fxx).*abs(-gxx./dzdx + d2logxdz2 - g_logw).^p;
    %pre_score = exp(-fxx).*abs(-gxx./dzdx + d2logxdz2 - g_logw);
    score = log(sum(pre_score)*(10/length(xx)))/p;
    
    figure(4)
    subplot(3,3,i)
    plot(z, pre_score)
    title([num2str(i), ' ', num2str(score)])
    set(gca, 'xlim', [-1,1])
end

zz = linspace(-1,1,1E4);
logw_zz = eval_log_measure(poly, zz);
B_zz = eval_basis(poly,zz);
figure
for i = 1:numel(maps)
    x_nodes = reference2domain(maps{i}, poly.nodes);
    dlogxdz = reference2domain_log_density(maps{i}, poly.nodes);
    logw = eval_log_measure(poly, poly.nodes);
    
    f = fun(x_nodes);
    coeff = poly.node2basis*exp(-0.5*(f(:)-dlogxdz(:)+logw(:)));
    %
    tmp = cumsum(B_zz.*coeff(:)',2);
    
    xx = reference2domain(maps{i}, zz);
    ff = fun(xx);
    dlogxxdz = reference2domain_log_density(maps{i}, zz);
    fe = exp(-0.5*(ff(:)-dlogxxdz(:)+logw_zz(:)));
    
    figure(5)
    subplot(3,3,i)
    semilogy(0:poly.order, sum((tmp-fe).^2,1))
    title(num2str(i))
    
    figure(6)
    subplot(3,3,i)
    plot(zz(:), tmp(:,[ceil(poly.order/2), poly.order+1]))
    hold on
    plot(zz(:), fe)
    title(num2str(i))
end

end