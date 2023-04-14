% test the interpolation

% first test the global one
nums = 2:11;
lags = cell(1, length(nums));

for j = 1:length(nums)
    lags{1, j} = Lagrangep(nums(j), 1);
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
hold on
lag = lags{1,end};
I   = eye(cardinal(lag));
lf  = zeros(size(xs));
for j = 1:cardinal(lag)
    fi  = eval(lag, I(:,j), xs);
    lf  = lf + abs(fi);
    plot(xs, fi)
    plot(lag.nodes, I(:,j), 'o')
end
plot(xs, lf, 'linewidth', 2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def = Lagrangep(5, 5);
figure
plot(eval_basis(def, def.nodes))

fn = figure;
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), figure(fn), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%

def = Lagrangep(10, 3);
figure
plot(eval_basis(def, def.nodes))

fn = figure;
xs = linspace(def.domain(1), def.domain(2), 1000);
bs = eval_basis(def, xs);
for i = 1:size(bs,2), figure(fn), plot(xs, bs(:,i)); set(gca,'ylim', [-1,1]); title(num2str(i)); pause(0.2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


def = Lagrangep(10, 5);
xs = linspace(def.domain(1), def.domain(2), 1000);
% test interpolation property
f = @(y) sin(y*5*pi)+1;
fi = eval(def, f(def.nodes), xs);
figure
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



def = Lagrangep(10, 5);
xs = linspace(def.domain(1), def.domain(2), 1000);
f = @(y) 1./log(2+y/2) + sin(y*5*pi);
fi = eval(def, f(def.nodes), xs);
figure
plot(xs, f(xs), xs, fi, def.nodes,f(def.nodes(:)), 'o')

%%%%
tic;
bs = eval_basis(def, xs);
for i = 1:10000
    f1 = bs*rand(cardinal(def),1);
end
toc
tic;
bs = sparse(eval_basis(def, xs));
for i = 1:10000
    f1 = bs*rand(cardinal(def),1);
end
toc
%
tic;
for i = 1:10000
    bs = eval_basis(def, xs);
    f1 = bs*rand(cardinal(def),1);
end
toc
tic;
for i = 1:10000
    f = eval(def, rand(cardinal(def),1), xs);
end
toc

%%%% test integration

err = zeros(5,6);

f = @(y) 1./log(y/2+2) + sin(y*5*pi);
a = integral(f, -1,1);
for i = 1:5
    for j = 1:6
        def = Lagrangep(2+j, 1+i);
        err(i,j) = abs(def.unweighed_int_W*f(def.nodes(:)) - a);
    end
end
figure
semilogy(err')
