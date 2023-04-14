
%%%% step 1: generate a 2D data set, scatter plot and histograms
d = 2; n = 1E5;
data = randn(d,d)*(random('beta',2,3,d,n).*randn(d,n));
figure
subplot(2,4,[1 2 5 6])
scatter(data(1,:), data(2,:))
xlabel('x_1')
ylabel('x_2')
subplot(2,4,[3 4])
histogram(data(1,:))
xlabel('x_1')
subplot(2,4,[7 8])
histogram(data(2,:))
xlabel('x_2')


%%%% step 2: build the empirical map and plot marginal pdfs and cdfs
%%%% the bound for the truncated Gaussian is [-2, 2], default is 4
EM = EmpiricalMap(data, 'gauss_bound', 2);

figure
subplot(2,2,1)
plot_pdf(EM,1);
title('pdf x_1')
subplot(2,2,3)
plot_cdf(EM,1);
title('cdf x_1')
%
subplot(2,2,2)
plot_pdf(EM,2);
title('pdf x_2')
subplot(2,2,4)
plot_cdf(EM,2);
title('cdf x_2')

%%%% step 3: transform the data using the cdf MAP, scatter plot and hist
u = eval_cdf(EM, data);
figure
subplot(2,4,[1 2 5 6])
scatter(u(1,:), u(2,:))
xlabel('u_1')
ylabel('u_2')
subplot(2,4,[3 4])
histogram(u(1,:))
xlabel('u_1')
subplot(2,4,[7 8])
histogram(u(2,:))
xlabel('u_2')

%%%% step 4: check the accuracy of the inverse map
y = eval_icdf(EM, u);
norm(y-data,'fro')

%%%% step 5: KDE plots of the original data and the transform
figure
subplot(1,2,1)
ksdensity(data', 'PlotFcn', 'contour')
title('data')
subplot(1,2,2)
ksdensity(u', 'PlotFcn', 'contour')
title('transformed')


%%%% step 6: transform the data using the Gauss MAP, scatter plot and hist
u = eval_map(EM, data);
figure
subplot(2,4,[1 2 5 6])
scatter(u(1,:), u(2,:))
xlabel('u_1')
ylabel('u_2')
subplot(2,4,[3 4])
histogram(u(1,:))
xlabel('u_1')
subplot(2,4,[7 8])
histogram(u(2,:))
xlabel('u_2')

%%%% step 7: check the accuracy of the inverse map
z = eval_imap(EM, u);
norm(z-data,'fro')

%%%% step 8: KDE plots of the original data and the transform
figure
subplot(1,2,1)
ksdensity(data', 'PlotFcn', 'contour')
title('data')
subplot(1,2,2)
ksdensity(z', 'PlotFcn', 'contour')
title('transformed')

