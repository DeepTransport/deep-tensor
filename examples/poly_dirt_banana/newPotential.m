function f = newPotential(func, beta, sirt, z)

[x, px] = sirt.eval_irt(z);
f = func(x).*beta - px;

end
