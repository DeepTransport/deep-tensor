function [core,interp_x,core_next] = build_basis_svd(oned, ...
    interp_xold, core_next, F, ...
    dir, int_method, loc_err_tol, max_rank)
% find the eigen function of the Schmidt operator
% cross case, nileft x niright kernels,
% dimension of the refinement block
nbleft  = size(F, 1);
nnodes  = size(F, 2);
nbright = size(F, 3);
m       = size(F, 4);
% dimension of the next core
rn1 = size(core_next, 1);
nn  = size(core_next, 2);
rn2 = size(core_next, 3);
%size(core_next)
%
if dir > 0
    % from left >k is integrated
    % T*T' give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F, nbleft, nnodes, nbright, []), [2,1,3,4]), nnodes*nbleft,  []); %(nright+enrich)*m
    rold    = nbleft;
    %rnext   = rn1;
else
    % <k is integrated
    % T'*T give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F, nbleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []); %(nbleft+enrich)*m
    rold    = nbright;
    %rnext   = rn2;
end
%
[B,A,r]     = TTFun.local_truncate(loc_err_tol, 1, max_rank, oned, F);
[ind,core,interp_atx] = point_selection(oned, int_method, B);
couple      = reshape(interp_atx*A, r, [], m);
interp_x    = TTFun.local_index(oned, dir, interp_xold, ind);
%
if dir > 0
    % from left >k is integrated
    % get the index for transpose A(:,:,j) for each j
    core = permute(reshape(core,nnodes,rold,r), [2,1,3]);
    %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m
    %block to the next block, first permute the block to m x r x rn1
    core_next = reshape(permute(couple(:,1:rn1,:), [3,1,2]), [], rn1) * reshape(core_next, rn1, []);
    %size(core_next)
    %[m, r, nn, rn2]
    core_next = permute(reshape(core_next, m, r, nn, rn2), [2, 3, 4, 1]);
else
    % <k is integrated
    % unfold the 3D tensor
    core    = permute(reshape(core,nnodes,rold,r), [3,1,2]);
    %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m
    %block to the next block, first permute the block to rn2 x r x m
    core_next = reshape(core_next, [], rn2) * reshape(permute(couple(:,1:rn2,:), [2,1,3]), rn2, []);
    core_next = reshape(core_next, rn1, nn, r, m);
end
end