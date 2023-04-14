function f = eval_ou_process_marginal(data, ind, x)

C = data.C(ind, ind);
f = 0.5*sum((C\x).*x,1);
z = 0.5*log(det(C)) + 0.5*length(ind)*log(2*pi);
f = z + f;

end