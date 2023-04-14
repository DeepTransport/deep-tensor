function f = newPotential(func, sirt, z)

[x, px] = sirt.eval_irt(z);
f = func(x) - px;

end
