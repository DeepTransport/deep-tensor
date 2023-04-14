function I = hyperbolic_cross(d, k, weights)

if nargin == 2
    weights = ones(1,d);
end
array = hyperbolic_cross_recur(d, k, weights);

I = MultiIndices(unique(array,'rows'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = hyperbolic_cross_recur(d, k, weights)

if d == 1
    I = reshape(0:k, [], 1);
else
    I_pre = hyperbolic_cross_recur(d-1,k,weights);
    n = size(I_pre,1);
    I = zeros(n*(k+1),d);
    %
    pos = 0;
    for j = 1:n
        tmp = prod(max(I_pre(j,:), 1)./weights(1:d-1));
        for i = 1:k+1
            if tmp*(max(i-1,1)/weights(d)) <= k
                pos = pos + 1;
                I(pos,1:d-1) = I_pre(j,:);
                I(pos,d) = i-1;
            end
        end
    end
    I = I(1:pos,:);
end

end
