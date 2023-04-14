function ess2n = ess_ratio(log_weight)

log_weight = log_weight - max(log_weight);
ess2n = ( sum(exp(log_weight)).^2/sum(exp(2*log_weight)) ) / length(log_weight(:));

end