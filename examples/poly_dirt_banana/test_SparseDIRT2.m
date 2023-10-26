clear all
%addpath(genpath('.'));

%% Generate the target = banana
d = 2;
xmin = -2;
xmax = 2;
nx = 100;

gridx = linspace(xmin,xmax,nx);
[gridxx,gridyy] = meshgrid(gridx,gridx);
X = [gridxx(:),gridyy(:)]';

sigma = 5;
a = 1;
model = Banana(sigma, a);
func = @(x,mask) model.eval_potential_dirt(x);

% plot
clf
fplot = reshape(func(X),nx,nx);
contourf(gridx,gridx,exp(-fplot))

%% 2-layers DIRT 

p = 30;
dom = BoundedDomain([xmin,xmax]);
base = ApproxBases(Legendre(p), dom, d);
%base = ApproxBases(Chebyshev1st(p), dom, d);

opt = SparseOption();
opt.tol = 1e-2;
opt.max_sample_size = 5e3;
opt.enrich_degree = 1;
opt.adaptation_rule = 'reducedmargin';
%
opt.init_sample_size = 2;
opt.enrich_sample_size = 0.5; 
%opt.weight_rule = 'none';
opt.fast = true;

beta = 0.2;
%sirt = SparseSIRT( @(x)func(x).*beta, base, s, 'defensive', 1E-15);
sirt = SparseSIRT( @(x)func(x).*beta, base, opt, 'defensive', 1E-5);

dom = BoundedDomain([0,1]);
base = ApproxBases(Legendre(p), dom, d);
%base = ApproxBases(Chebyshev1st(p), dom, d);
%newSirt = SparseSIRT( @(z)newPotential(func,1,sirt,z), base, s, 'defensive', 1E-10);


newSirt = SparseSIRT( @(z)newPotential(func,1,sirt,z), base, opt, 'defensive', 1E-5);


%%% plot

clf

%% Plot function, layer 1
subplot(2,2,1)
f1 = func(X)*beta;
fplot = reshape(f1,nx,nx);
contourf(gridx,gridx,exp(-fplot))
title('Target (unormalized) density')
colorbar

subplot(2,2,2)
f2 = sirt.eval_potential(X);
fplot = reshape(f2,nx,nx);
contourf(gridx,gridx,exp(-fplot))
title('SparseSIRT')
colorbar

subplot(2,2,3)
fplot = reshape(f2-f1,nx,nx);
contourf(gridx,gridx,fplot)
title('Log(targetPDF/SIRT)')
colorbar

subplot(2,2,4)
gridz = linspace(0,1,nx);
[gridz1,gridz2] = meshgrid(gridz,gridz);
Z = [gridz1(:),gridz2(:)]';
[RZ, PZ] = sirt.eval_irt(Z);

pullbackTar = PZ-func(RZ)*beta;

fplot = reshape( pullbackTar ,nx,nx);
contourf(gridz,gridz,fplot)
title('Log of the pullback target')
colorbar


%% Plot grid, layer 1
figure
nz = 100;
gridz = linspace(0,1,nz);
[gridz1,gridz2] = meshgrid(gridz,gridz);
Z = [gridz1(:),gridz2(:)]';

subplot(1,2,1)
plot(Z(1,:),Z(2,:),'.')
title('Grid on the reference domain')

subplot(1,2,2)
TZ = sirt.eval_irt(Z);
plot(TZ(1,:),TZ(2,:),'.')
title('Image grid on the domain')


%% Plot function, layer 2

nz = 100;
gridz = linspace(0,1,nz);
[gridz1,gridz2] = meshgrid(gridz,gridz);
Z = [gridz1(:),gridz2(:)]';

figure

subplot(1,2,1)
f1 = newPotential(func,1,sirt,Z);
fplot = reshape(f1,nz,nz);
contourf(gridz,gridz,exp(-fplot))
title('Target (unormalized) density')
colorbar

subplot(1,2,2)
f2 = newSirt.eval_potential(Z);
fplot = reshape(f2,nz,nz);
contourf(gridx,gridx,exp(-fplot))
title('SparseSIRT')
colorbar


%% Plot grid, layer 2
figure

subplot(1,3,1)
plot(Z(1,:),Z(2,:),'.')
title('Grid on the reference domain')

subplot(1,3,2)
T1Z = newSirt.eval_irt(Z);
plot(T1Z(1,:),T1Z(2,:),'.')
title('Image of the grid by NewSIRT')

subplot(1,3,3)
TZ = sirt.eval_irt(T1Z);
plot(TZ(1,:),TZ(2,:),'.')
title('Image of the grid in the domain')
