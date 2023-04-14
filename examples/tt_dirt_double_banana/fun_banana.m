function [mllkd, mlp, gmllkd, gmlp] = fun_banana(u, data, sigma, beta)

F   = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
mllkd = sum((F-data).^2,1)*beta/(2*sigma^2);
mlp = 0.5*sum(u.^2,1);

%lpt = -sum((F-data).^2,1)*beta/(2*sigma^2) - sum(u.^2,1)*beta/2;%
%p   = exp(lpt);

tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
    ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );

gmllkd = (F*numel(data) - sum(data(:))).*tmp/sigma^2;
gmlp = u;


end

%{
sig = 0.3;
%dat = [3.940392199905546; 4.403271869259551]; % 2 bananas
dat = [3; 5]; % 2 bananas
beta = 1;

n  = 50;
xs = linspace(-4, 4, n);
ys = linspace(-4, 4, n);
[xx,yy] = meshgrid(xs, ys);
xts = [xx(:), yy(:)]';


[mllkd, mlp, gmllkd, gmlp] = fun_banana(xts, dat, sig, beta);

mllkd_p = gmllkd;
mllkd_m = gmllkd;
mlp_p = gmlp;
mlp_m = gmlp;

tol = 1E-5;
xts_p = xts;
xts_p(1,:) = xts_p(1,:)+tol;
[mllkd_p(1,:), mlp_p(1,:)] = fun_banana(xts_p, dat, sig, beta);
xts_p = xts;
xts_p(2,:) = xts_p(2,:)+tol;
[mllkd_p(2,:), mlp_p(2,:)] = fun_banana(xts_p, dat, sig, beta);

xts_m = xts;
xts_m(1,:) = xts_m(1,:)-tol;
[mllkd_m(1,:), mlp_m(1,:)] = fun_banana(xts_m, dat, sig, beta);
xts_m = xts;
xts_m(2,:) = xts_m(2,:)-tol;
[mllkd_m(2,:), mlp_m(2,:)] = fun_banana(xts_m, dat, sig, beta);

gmllkd_d = (mllkd_p - mllkd_m)/(2*tol);
gmlp_d = (mlp_p - mlp_m)/(2*tol);

norm(gmllkd(:) - gmllkd_d(:))
norm(gmlp(:) - gmlp_d(:))


%}