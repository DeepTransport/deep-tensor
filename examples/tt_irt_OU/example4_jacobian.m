
z = rand(d, 1E2);
for i = 1:size(irts,1)
    for j = 1:min(2,size(irts,2))
        irts{i,j} = set_defensive(irts{i,j}, 5E-2);
        debug_jac(irts{i,j}, z, 1);
        debug_jac(irts{i,j}, z, -1);
    end
end

