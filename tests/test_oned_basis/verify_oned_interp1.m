% test the interpolation

% first test the global one
nums = 4:3:31;
lags = cell(1, length(nums));

for j = 1:length(nums)
    lags{1,j} = Lagrange1(nums(j));
end

%%%%%%%% case 1

test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) sin(y*5*pi))
test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) sin(y*5*pi))

test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) sqrt(y+1));
test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) sqrt(y+1));

test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) 1./sqrt(y+1));
test_interp_conv(lags(1,[5, 7, 8, 10]), @(y)  1./sqrt(y+1));

test_interp_conv(lags(1,[1, 2, 4, 5]), @(y) log(y+1));
test_interp_conv(lags(1,[5, 7, 8, 10]), @(y) log(y+1));

%%%%%%%% case 2

xs  = linspace(-1,1,1000)';
figure
lag = lags{1,end};
I   = eye(cardinal(lag));
lf  = zeros(size(xs));
for j = 1:cardinal(lag)
    fi  = eval_radon(lag, I(:,j), xs);
    lf  = lf + abs(fi);
    plot(xs, fi)
    hold on
    plot(lag.nodes, I(:,j), 'o')
end
plot(xs, lf, 'linewidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
def = Lagrange1(15);
plot(eval_basis(def, def.nodes))

fn = figure;
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), figure(fn), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%
figure
def = Lagrange1(21);
plot(eval_basis(def, def.nodes))

fn = figure;
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), figure(fn), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def = Lagrange1(50);

xs = linspace(def.domain(1), def.domain(2), 1000);
% test interpolation property
f = @(y) sin(y*5*pi)+1;
fi = eval_radon(def, f(def.nodes), xs);
figure
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def = Lagrange1(50);
xs = linspace(def.domain(1), def.domain(2), 1000);
f = @(y) 1./log(2+y/2) + sin(y*5*pi);
fi = eval_radon(def, f(def.nodes(:)), xs);
figure
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%
f = zeros(1000, 10000);
A = rand(cardinal(def),10000);
tic;
bs = eval_basis(def, xs);
for i = 1:10000
    f(:,i) = bs*A(:,i);
end
toc
tic;
bs = eval_basis(def, xs);f = bs*A;toc
tic;
for i = 1:10000
    bs = eval_basis(def, xs);
    f(:,i) = bs*A(:,i);
end
toc
tic;
for i = 1:10000
    f(:,i) = eval_radon(def, A(:,i), xs);
end
toc
tic;f = eval_radon(def, A, xs);toc

%%%% test integration

err = zeros(1,6);

f = @(y) 1./log(2+y/2) + sin(y*5*pi);
a = integral(f,-1,1);
for j = 1:6
    def = Lagrange1(5*j);
    err(j) = abs(def.unweighed_int_W*f(def.nodes(:)) - a);
end
semilogy(err')
