function f = eval_ou_process(data, x)

f = 0.5*sum((data.B*x).^2,1);

end