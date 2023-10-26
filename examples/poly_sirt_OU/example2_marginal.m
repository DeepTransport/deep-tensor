% Does not pass the test

% should test ind = 1, ind = 1:(d-1) for > 0
% should test ind = d, ind = 2:d for < 0
% sample
z = rand(d, 1E4);
for i = 1:length(bases)
    figure;
    % test 1
    ind  = 1:1;
    if irts{i}.int_dir ~= 1, irts{i} = marginalise(irts{i}, 1); end
    tic;[r,p] = eval_irt(irts{i}, z(ind,:));toc
    fx = eval_pdf(irts{i}, r);
    tic;z0 = eval_rt(irts{i}, r);toc
    disp(' ')
    disp(['transform eror: ' num2str(norm(z(ind,:) - z0))])
    disp(['density eror: ' num2str(norm(exp(-p) - fx))])
    fe = eval_potential_marginal(data, ind, r);
    disp(['approx eror: ' num2str(norm(exp(-p) - exp(-fe)))])
    disp(' ')
    %
    subplot(2,3,1);plot(abs(fe - p)/max(abs(fe)), '.');
    subplot(2,3,2);plot(fe , p, '.');
    title('actual poential function value vs fft')
    subplot(2,3,3);plot(data.C(ind, ind) - cov(r'))
    
    % test 2
    ind  = 1:2;
    if irts{i}.int_dir ~= -1, irts{i} = marginalise(irts{i}, -1); end
    tic;[r,p] = eval_irt(irts{i}, z(ind,:));toc
    fx = eval_pdf(irts{i}, r);
    tic;z0 = eval_rt(irts{i}, r);toc
    disp(' ')
    disp(['transform eror: ' num2str(norm(z(ind,:) - z0))])
    disp(['density eror: ' num2str(norm(exp(-p) - fx))])
    %
    fe = eval_potential_marginal(data, ind, r);
    disp(['approx eror: ' num2str(norm(exp(-p) - exp(-fe)))])
    disp(' ')
    subplot(2,3,4);plot(abs(fe - p)/max(abs(fe)), '.');
    subplot(2,3,5);plot(fe , p, '.');
    title('actual potential function value vs fft')
    subplot(2,3,6);plot(data.C(ind, ind) - cov(r'))
end



