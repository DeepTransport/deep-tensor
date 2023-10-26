% sample
z = rand(d, 1E4);
for i = 1:size(irts,1)
    for j = 1:size(irts,2)
        figure
        tic;[r,f] = eval_irt(irts{i,j}, z);toc
        tic;z0 = eval_rt(irts{i,j}, r);toc
        disp(' ')
        disp(['transform eror: ' num2str(norm(z-z0, 'fro'))])
        disp(['potential eror: ' num2str(norm(func(r) - f))])
        disp(['pdf eror: ' num2str(norm(exp(-func(r)) - exp(-f)))])
        disp(['cov eror: ' num2str(norm(data.C - cov(r'))/norm(data.C))])
        disp(' ')
        %
        subplot(2,2,1); plot(abs(func(r) - f), '.'); title('actual potential function vs fft')
        subplot(2,2,2); plot(func(r), f, '.');
        subplot(2,2,3); plot(data.C - cov(r')); title('actual covariance vs sample covariance')
    end
end

