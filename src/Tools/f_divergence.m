function [dkl, dh2, dtv] = f_divergence(lp_ref, lp)

n = length(lp_ref);
m = size(lp, 1);

lp_ref = repmat(lp_ref, m, 1);

t = lp - lp_ref;
%
logratio = logavrexp(t); % ratio of normalising const p over p_ref
% 
dkl = ( sum(-t, 2)/n + logratio(:) )';
% 
dh2 = 1 - exp( logavrexp(0.5*t) - 0.5*logratio )';
%
dtv = 0.5*sum(abs(exp(t + logratio(:)) - 1), 2)/n;

end

function a = logavrexp(t)

[m, n] = size(t);
a = zeros(m, 1);

for i = 1:m
    mt   = max(t(i,:));
    a(i) = mt + log(sum(exp(t(i,:)-mt))) - log(n);
end

end