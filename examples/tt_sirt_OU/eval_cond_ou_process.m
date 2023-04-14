function f = eval_cond_ou_process(data, x1, x2, dir)

%order of the samples (x1, x2)
x = [x1; x2];
f = 0.5*sum((data.B*x).^2,1) + log(data.norm);

if dir > 0
    %from the last dimension, marginalise to the first
    %conditional (x2 | x1)
    ind = 1:size(x1,1);
    fm  = eval_ou_process_marginal(data, ind, x1);
else
    %from the first dimension, marginalise to the last
    %conditional (x1 | x2)
    ind = size(x1, 1) + ( 1:size(x2,1) );
    fm  = eval_ou_process_marginal(data, ind, x2);
end

f = f - fm;

end