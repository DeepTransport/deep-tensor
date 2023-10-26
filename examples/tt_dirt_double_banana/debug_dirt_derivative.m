% need to run example dirt first to obtain airt and eirt

n = 1E2;
r = random(airt.ref, d, n);

[xa, fa, ga] = eval_irt(airt, r);
[xe, fe, ge] = eval_irt(eirt, r);

tol = 1E-5;

norm(eval_potential(airt, xa) - fa)
norm(eval_potential(eirt, xe) - fe)

fp = zeros(d,n);
fm = zeros(d,n);
for i = 1:d
    rp = r;
    rp(i,:) = rp(i,:) + tol;
    [~,fp(i,:)] = eval_irt(airt, rp);
    rm = r;
    rm(i,:) = rm(i,:) - tol;
    [~,fm(i,:)] = eval_irt(airt, rm);
end
gad = (fp-fm)/(2*tol);
norm(ga(:) - gad(:))/norm(ga(:))

fp = zeros(d,n);
fm = zeros(d,n);
for i = 1:d
    rp = r;
    rp(i,:) = rp(i,:) + tol;
    [~,fp(i,:)] = eval_irt(eirt, rp);
    rm = r;
    rm(i,:) = rm(i,:) - tol;
    [~,fm(i,:)] = eval_irt(eirt, rm);
end
ged = (fp-fm)/(2*tol);
norm(ge(:) - ged(:))/norm(ge(:))

