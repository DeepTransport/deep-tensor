
% setup the OU process
d = 20;
a = 0.5;
%p = randperm(d,d);
%pt(p) = 1:length(p);


% test m = 1, m = d-1, m = 8
% sample
z = rand(d, 1E4);
for i = 1:size(irts,1)
    for j = 1:min(2,size(irts,2))
        m = 8;
        ind1 = 1:m;
        ind2 = (m+1):d;
        xc1 = debug_x(ind1,:);
        xc2 = debug_x(ind2,:);
        % conditional sample
        zc1 = z(ind1, :);
        zc2 = z(ind2, :);
        %
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[rc2,f] = eval_cirt(irts{i,j}, xc1, zc2);toc
        %
        figure;
        subplot(2,2,1);plot(abs(exp(-eval_potential_conditional(data, xc1, rc2, 1)) - exp(-f))); title('actual function value vs fft')
        subplot(2,2,2);plot(data.C(ind2, ind2) - cov(rc2')); title('actual covariance vs sample covariance')
        %
        disp(' ')
        disp(['approx eror: ' num2str(norm(exp(-eval_potential_conditional(data, xc1, rc2, 1)) - exp(-f)))])
        disp(['cov eror: ' num2str(norm(data.C(ind2, ind2) - cov(rc2'))/norm(data.C))])
        disp(' ')
        
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[rc1,f] = eval_cirt(irts{i,j}, xc2, zc1);toc
        %
        subplot(2,2,3);plot(abs(exp(-eval_potential_conditional(data, rc1, xc2, -1)) - exp(-f))); title('actual function value vs fft')
        subplot(2,2,4);plot(data.C(ind1, ind1) - cov(rc1')); title('actual covariance vs sample covariance')
        disp(' ')
        disp(['approx eror: ' num2str(norm(exp(-eval_potential_conditional(data, rc1, xc2, -1)) - exp(-f)))])
        disp(['cov eror: ' num2str(norm(data.C(ind1, ind1) - cov(rc1'))/norm(data.C))])
        disp(' ')
    end
end

for i = 1:size(irts,1)
    for j = 1:min(2,size(irts,2))
        irts{i,j} = set_defensive(irts{i,j}, 1E-2);
        m = 8;
        indx = 1:m;
        indy = (m+1):d;
        x = debug_x(indx,1);
        my = data.C(indy,indx)*(data.C(indx,indx)\x); 
        Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
        % conditional sample
        zy = z(indy, :);
        %
        if irts{i,j}.int_dir ~= 1, irts{i,j} = marginalise(irts{i,j}, 1); end
        tic;[y,f] = eval_cirt(irts{i,j}, x, zy);toc
        tic;
        zx = eval_rt (irts{i,j}, x);
        fx = eval_potential(irts{i,j}, x);
        [r,fxy] = eval_irt(irts{i,j}, [repmat(zx, 1, size(zy,2)); zy]);
        y2 = r(indy,:);
        f2 = fxy - fx;
        toc
        disp([norm(y2-y, 'fro'), norm(f2-f)])
        
        %
        figure;
        subplot(2,3,1);plot(abs(exp(-eval_potential_conditional(data, repmat(x,1,size(zy,2)), y, 1)) - exp(-f))); title('actual function value vs fft')
        subplot(2,3,2);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
        subplot(2,3,3);plot(my - mean(y,2)); title('actual mean vs sample mean')
        
        disp(' ')
        disp(['approx eror: ' num2str(norm(exp(-eval_potential_conditional(data, repmat(x,1,size(zy,2)), y, 1)) - exp(-f)))])
        disp(['cov eror: ' num2str(norm(Cy - cov(y'))/norm(data.C))])
        disp(' ')
        
        m = 8;
        indy = 1:m;
        indx = (m+1):d;
        x = debug_x(indx,1);
        my = data.C(indy,indx)*(data.C(indx,indx)\x); 
        Cy = data.C(indy,indy) - data.C(indy,indx)*(data.C(indx,indx)\data.C(indx,indy));
        % conditional sample
        zy = z(indy, :);
        %
        if irts{i,j}.int_dir ~= -1, irts{i,j} = marginalise(irts{i,j}, -1); end
        tic;[y,f] = eval_cirt(irts{i,j}, x, zy);toc
        %
        subplot(2,3,4);plot(abs(exp(-eval_potential_conditional(data, y, repmat(x,1,size(zy,2)), -1)) - exp(-f))); title('actual function value vs fft')
        subplot(2,3,5);plot(Cy - cov(y')); title('actual covariance vs sample covariance')
        subplot(2,3,6);plot(my - mean(y,2)); title('actual mean vs sample mean')
        
        disp(' ')
        disp(['approx eror: ' num2str(norm(exp(-eval_potential_conditional(data, y, repmat(x,1,size(zy,2)), -1)) - exp(-f)))])
        disp(['cov eror: ' num2str(norm(Cy - cov(y'))/norm(data.C))])
        disp(' ')
    end
end
