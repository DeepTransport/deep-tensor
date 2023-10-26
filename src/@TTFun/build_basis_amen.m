function [core, interp_x, res_w, res_x, core_next] = build_basis_amen(oned, ...
    interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr, ...
    dir, int_method, loc_err_tol, max_rank, kick_rank)
% find the eigen function of the Schmidt operator
% cross case, nileft x niright kernels,
%
nbleft  = size(F, 1);
nbright = size(F, 3);
nrleft  = size(Fr,1);
nrright = size(Fr,3);
nnodes  = cardinal(oned);
m       = size(F, 4);
% dimension of the next core
rn1 = size(core_next, 1);
nn  = size(core_next, 2);
rn2 = size(core_next, 3);
%
if dir > 0
    % from left >k is integrated
    % T*T' give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,1,3,4]), nnodes*nbleft, []);
    Fu      = reshape(permute( reshape(Fu, nbleft, nnodes, nrright, []), [2,1,3,4]), nnodes*nbleft, []);
    rold    = nbleft;
    rnext   = rn1;
else
    % <k is integrated
    % T'*T give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
    Fu      = reshape(permute( reshape(Fu, nrleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
    rold    = nbright;
    rnext   = rn2;
end
%
[B,A,r] = TTFun.local_truncate(loc_err_tol, 1, max_rank, oned, F);
%
% compute residual
% BSV' = FT is arranged as nodes, rold (fixed), rnew (to be updated)
% rnew dimension needs to be projected onto res basis from the update
% direction, rold dimension needs to be projected onto res basis from
% the fixed dimension
% for the right projection
if dir > 0
    % A is r x rnext x m, multiply the rnext dim with res_w_r
    tmp_r   = reshape(permute(reshape(A,r,rn1,m), [1,3,2]), [], rn1)*res_w_r;
    tmp_r   = reshape(permute(reshape(tmp_r,r,m,[]),[1,3,2]),r,[]);
    Fu      = Fu - B*tmp_r;
    % for the left projection
    tmp_l   = res_w_l*reshape(permute(reshape(B,nnodes,rold,r), [2,1,3]), rold, []);
    tmp_l   = reshape(tmp_l, nrleft*nnodes, r);
    % align Fr as rold (nrleft), nodes, rnew (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,[]) - reshape(tmp_l*tmp_r,nrleft,nnodes,nrright,[]);
    Fr      = reshape(permute(Fr,[2,1,3,4]), nnodes*nrleft, []);
    rrold   = nrleft;
else
    % A is r x rnext x m, multiply the rnext dim with res_w_l
    tmp_lt  = res_w_l*reshape(permute(reshape(A,r,rn2,m), [2,1,3]), rn2, []);
    % r x nr_l x m
    tmp_lt  = reshape(permute(reshape(tmp_lt,[],r,m),[2,1,3]),r,[]);
    % Fu is aligned with F
    Fu      = Fu - B*tmp_lt;
    % for the right projection
    tmp_r   = reshape(permute(reshape(B,nnodes,rold,r), [3,1,2]), [], rold)*res_w_r; % r x n x nr_r
    tmp_r   = reshape(tmp_r, r, nnodes*nrright);
    % align Fr as rnew (nrleft), nodes, rold (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,m) - permute(reshape(tmp_lt'*tmp_r,nrleft,m,nnodes,nrright), [1,3,4,2]);
    Fr      = reshape(permute(Fr,[2,3,1,4]), nnodes*nrright, []);
    rrold   = nrright;
end
% enrich basis
if isa(oned, 'Piecewise')
    % T = kron(speye(rold), oned.mass_L') * [B, Fu];
    T = reshape(oned.mass_R*reshape([B,Fu],cardinal(oned),[]),size(B,1),[]);
    [Q,R] = qr(T,0);
    if size(T,1) < size(T,2)
        Q = Q(:,1:size(T,1));
    end
    % interpolation basis
    % B = full(kron(speye(rold),oned.mass_L') \ Q);
    B = reshape(oned.mass_R\reshape(Q,cardinal(oned),[]),size(Q,1),[]);
else
    T = [B,Fu];
    [B,R] = qr(T,0);
    if size(T,1) < size(T,2)
        B = B(:,1:size(T,1));
    end
end
rnew        = size(B,2);
%
[indf,core,interp_atx] = point_selection(oned, int_method, B);
couple      = reshape(interp_atx*(R(1:rnew,1:r)*A), rnew, rnext, m);

%{
[B,A,r] = TTFun.local_truncate(loc_err_tol, 1, max_rank, oned, F);
if dir > 0
    % A is r x rnext x m, multiply the rnext dim with res_w_r
    tmp_r   = reshape(permute(reshape(A,r,rn1,m), [1,3,2]), [], rn1)*res_w_r;
    tmp_r   = reshape(permute(reshape(tmp_r,r,m,[]),[1,3,2]),r,[]);
    % for the left projection
    tmp_l   = res_w_l*reshape(permute(reshape(B,nnodes,rold,r), [2,1,3]), rold, []);
    tmp_l   = reshape(tmp_l, nrleft*nnodes, r);
    % align Fr as rold (nrleft), nodes, rnew (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,[]) - reshape(tmp_l*tmp_r,nrleft,nnodes,nrright,[]);
    Fr      = reshape(permute(Fr,[2,1,3,4]), nnodes*nrleft, []);
    rrold   = nrleft;
else
    % A is r x rnext x m, multiply the rnext dim with res_w_l
    tmp_lt  = res_w_l*reshape(permute(reshape(A,r,rn2,m), [2,1,3]), rn2, []);
    % r x nr_l x m
    tmp_lt  = reshape(permute(reshape(tmp_lt,[],r,m),[2,1,3]),r,[]);
    % for the right projection
    tmp_r   = reshape(permute(reshape(B,nnodes,rold,r), [3,1,2]), [], rold)*res_w_r; % r x n x nr_r
    tmp_r   = reshape(tmp_r, r, nnodes*nrright);
    % align Fr as rnew (nrleft), nodes, rold (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,m) - permute(reshape(tmp_lt'*tmp_r,nrleft,m,nnodes,nrright), [1,3,4,2]);
    Fr      = reshape(permute(Fr,[2,3,1,4]), nnodes*nrright, []);
    rrold   = nrright;
end
[B,A,r]     = TTFun.local_truncate(loc_err_tol, 1, max_rank, oned, [F,Fu]);
rnew        = size(B,2);
%
[indf,core,interp_atx] = point_selection(oned, int_method, B);
couple      = reshape(interp_atx*A, r, [], m);
%}
%
interp_x    = TTFun.local_index(oned, dir, interp_xold, indf);
%
Qr      = TTFun.local_truncate(loc_err_tol*eps, 1, kick_rank, oned, Fr);
% rold is the leading dimension in the point selection
indr    = point_selection(oned, int_method, Qr);
res_x   = TTFun.local_index(oned, dir, res_xold, indr);
%
if dir > 0
    % from left >k is integrated
    % get the index for transpose A(:,:,j) for each j
    core = permute(reshape(core,nnodes,rold,rnew), [2,1,3]);
    %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m
    %block to the next block, first permute the block to m x r x rn1
    core_next = reshape(permute(couple(:,1:rnext,:), [3,1,2]), [], rnext) * reshape(core_next, rnext, []);
    core_next = permute(reshape(core_next, m, rnew, nn, rn2), [2, 3, 4, 1]);
    %
    tmp     = reshape(permute(reshape(res_w_l*reshape(core, rold, nnodes*rnew),rrold,nnodes,rnew),[2,1,3]),[],rnew);
    res_w   = tmp(indr,:);
else
    % <k is integrated
    % unfold the 3D tensor
    core    = permute(reshape(core,nnodes,rold,rnew), [3,1,2]);
    %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m
    %block to the next block, first permute the block to rn2 x r x m
    core_next = reshape(core_next, [], rnext) * reshape(permute(couple(:,1:rnext,:), [2,1,3]), rnext, []);
    core_next = reshape(core_next, rn1, nn, rnew, m);
    %
    tmp     = reshape(reshape(core, [], rold)*res_w_r, rnew,[]);
    res_w   = tmp(:,indr);
end
end