function rerr = local_error(core, f)
% Compute local error for TT block
%
diff = core(:) - f(:);
rerr = max(abs(diff)) / max(abs(f(:)));
end