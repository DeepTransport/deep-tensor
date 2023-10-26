% Wrapper from my to FTT interface
function [mllkds, mlps] = diffusion_lpfun_ftt(theta, phil, bound, W2, Mass_summed, Q_obs, sigma_n, masks)
mllkds = ((diffusion_QoI(theta', phil, bound, W2, Mass_summed) - Q_obs).^2) / (2*sigma_n);

if (nargin<8) || (isempty(masks)) || (all(masks>0))
    if (nargin>7) && (~isempty(masks))
        mllkds = mllkds(:, masks);
    end
    mllkds = sum(mllkds, 2);
end

mllkds = mllkds';

% mlps = zeros(1, size(theta,2));
mlps = sum(theta.^2,1)/2;
end
