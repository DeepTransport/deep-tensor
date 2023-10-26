
z = rand(d, 123);
for i = 1:length(bases)
    irts{i} = set_defensive(irts{i}, 5E-6);
    debug_jac(irts{i}, z, [3 1 2]);
    debug_jac(irts{i}, z, [1 3 2]);
end

