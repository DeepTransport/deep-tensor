function debug_jac(irt, z, dir)

if irt.int_dir ~= dir, irt = marginalise(irt, dir); end
%
tic;ur = eval_irt(irt, z);toc
tic;zr = eval_rt (irt, ur);toc
%
tol = 1E-6;
tic;
[m,n] = size(ur);
ut = reshape(repmat(ur, m, 1), m, []) + repmat(eye(m)*tol, 1, n);
zt = eval_rt(irt, ut);
%
um = reshape(repmat(ur, m, 1), m, []) - repmat(eye(m)*tol, 1, n);
zm = eval_rt(irt, um);
%
Jd = ( zt - zm ) / (tol*2);
toc
%
tic; J = eval_rt_jac(irt, ur, zr); toc
%
disp(norm(J(:)-Jd(:)))

plot(J(:), Jd(:))

end

