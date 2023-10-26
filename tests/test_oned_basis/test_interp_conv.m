function test_interp_conv(lags, f)

xs  = linspace(-1,1,1000)';
fs  = f(xs);

figure
plot(xs, fs, 'linewidth', 2)
hold on

for j = 1:size(lags, 2)
    fx  = f(lags{j}.nodes);
    fi  = eval_radon(lags{j}, fx(:), xs);
    plot(xs, fi)
    plot(lags{j}.nodes, fx(:), 'o')
end


end