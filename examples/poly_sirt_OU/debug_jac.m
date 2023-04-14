function debug_jac(irt, z, dir)

irt = marginalise(irt, dir);
%
tic;ur = eval_irt_reference(irt, z);toc
tic;zr = eval_rt_reference (irt, ur);toc
%
tol = 1E-5;
tic;
[m,n] = size(ur);
ut = reshape(repmat(ur, m, 1), m, []) + repmat(eye(m)*tol, 1, n);
zt = eval_rt_reference(irt, ut);
%
um = reshape(repmat(ur, m, 1), m, []) - repmat(eye(m)*tol, 1, n);
zm = eval_rt_reference(irt, um);
%
Jd = ( zt - zm ) / (tol*2);
toc
%
tic; J = eval_rt_jac_reference(irt, ur, zr); toc
%
disp(norm(J(:)-Jd(:)))

plot(J(:), Jd(:))

end

