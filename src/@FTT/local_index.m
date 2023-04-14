function interp_x = local_index(oned, direction, interp_xold, ind)
% Build the local nested index set
%
if isempty(interp_xold)
    % interpolation points
    interp_x = reshape( oned.nodes(ind), 1, [] );
else
    % the basis coefficients is now an unfolded tensor, npoly x rold x rnew
    % the candidate index set ( x_old, nodes ), nodes is the leading dimension
    nnodes  = length(oned.nodes);
    rold    = size(interp_xold, 2);
    ipair   = [ reshape( repmat(1:rold, nnodes, 1), 1, nnodes*rold); repmat(1:nnodes, 1, rold) ];
    iselect = ipair(:, ind);
    if direction > 0
        interp_x = [ interp_xold(:, iselect(1,:)); reshape( oned.nodes(iselect(2,:)), 1, length(ind)) ];
    else
        interp_x = [ reshape( oned.nodes(iselect(2,:)), 1, length(ind)); interp_xold(:, iselect(1,:)) ];
    end
end
end