
xmin = -5;
xmax = 5;
nx = 100;

gridx = linspace(xmin,xmax,nx);
[gridxx,gridxy] = meshgrid(gridx,gridx);
X = [gridxx(:),gridxy(:)]';
%
nz = nx;
gridz = linspace(-1,1,nz);
[gridzx,gridzy] = meshgrid(gridz,gridz);
Z = [gridzx(:),gridzy(:)]';
%
nu = nx;
gridu = linspace(1E-2,1-1E-2,nu);
[gridux,griduy] = meshgrid(gridu,gridu);
U = [gridux(:),griduy(:)]';


%%
sigma = 2;
a = 1;

beta = 0.05;
% plot
figure
subplot(2,2,1)
fplot = reshape(banana(X,sigma,a),nx,nx);
contourf(gridx,gridx,exp(-fplot))
subplot(2,2,2)
fplot = reshape(banana(X,sigma/beta,a),nx,nx);
contourf(gridx,gridx,exp(-fplot))

d = 2;
p = 50;
c1 = Chebyshev1st(p);
c2 = Chebyshev2nd(p);
le = Legendre(p);
%

opt = SparseOption();
opt.tol = 1e-2;
opt.max_sample_size = 5e3;
opt.enrich_degree = 2;
opt.adaptation_rule = 'margin';

%{
bases = {};
for i = 1:8
    bases{i} = ApproxBases(c1, AlgebraicMapping(i), d);
end
scores = estimate_scaling(@(x)banana(x,sigma/beta,a), bases, 2, 1E4);
%

for s = [4,6,8,10]
    dom = AlgebraicMapping(s);
    base = ApproxBases(c1, dom, d);
    
    y = SIRT.potential_to_sqrt_density(base, @(x) banana(x,sigma,a), Z);
    subplot(2,2,3)
    ymapped = reshape(y,nz,nz);
    contourf(gridz,gridz,exp(-ymapped))
    
    y = SIRT.potential_to_sqrt_density(base, @(x) banana(x,sigma/beta,a), Z);
    subplot(2,2,4)
    ymapped = reshape(y,nz,nz);
    contourf(gridz,gridz,exp(-ymapped))
    
    sirt = SparseSIRT( @(x)banana(x, sigma/beta,a), base, opt, 'defensive', 1E-6);
end
%}

dom = AlgebraicMapping(6);
base = ApproxBases(le, dom, d); % c1

y = SIRT.potential_to_sqrt_density(base, @(x) banana(x,sigma,a), Z);
subplot(2,2,3)
ymapped = reshape(y,nz,nz);
contourf(gridz,gridz,exp(-ymapped))

y = SIRT.potential_to_sqrt_density(base, @(x) banana(x,sigma/beta,a), Z);
subplot(2,2,4)
ymapped = reshape(y,nz,nz);
contourf(gridz,gridz,exp(-ymapped))

sirt = SparseSIRT( @(x)banana(x, sigma/beta,a), base, opt, 'defensive', 1E-6);

%

newdom = BoundedDomain([0,1]);
newbase = ApproxBases(le, newdom, d);
newSirt = SparseSIRT( @(z)newPotential(@(x)banana(x,sigma,a),sirt,z), newbase, opt, 'defensive', 1E-5);

%%% plot

figure

subplot(2,2,1)
fplot = reshape(banana(X,sigma/beta,a),nx,nx);
contourf(gridx,gridx,exp(-fplot))
title('Target (unormalized) density')

subplot(2,2,2)
f1 = sirt.eval_potential(X);
fplot = reshape(f1,nx,nx);
contourf(gridx,gridx,exp(-fplot))
title('1st SparseSIRT')

subplot(2,2,3)
[RX, PX] = sirt.eval_irt(U);
pullbackTar = banana(RX,sigma,a)-PX;
fplot = reshape( pullbackTar ,nx,nx);
contourf(gridz,gridz,exp(-fplot))
title('the pullback target')

subplot(2,2,4)
f2 = newSirt.eval_potential(U);
fplot = reshape(f2,nx,nx);
contourf(gridx,gridx,exp(-fplot))
title('2nd SparseSIRT')


% Plot grid

figure 

subplot(2,2,1)
plot(U(1,:),U(2,:),'.')
title('Grid on the reference domain')

subplot(2,2,2)
TZ = sirt.eval_irt(U);
plot(TZ(1,:),TZ(2,:),'.')
hold on
fplot = reshape(banana(X,sigma/beta,a),nx,nx);
contour(gridx,gridx,exp(-fplot))
title('Image grid on the domain, 1st layer')


subplot(2,2,3)
T1Z = newSirt.eval_irt(U);
plot(T1Z(1,:),T1Z(2,:),'.')
title('Image of the grid by NewSIRT')

subplot(2,2,4)
TZ = sirt.eval_irt(T1Z);
plot(TZ(1,:),TZ(2,:),'.')
hold on
fplot = reshape(banana(X,sigma,a),nx,nx);
contour(gridx,gridx,exp(-fplot))
hold on
title('Image of the grid in the domain, 2nd layer')
