function I = total_degree(d, k)

I = MultiIndices(zeros(1,d));
for i = 1:k
    Iadd = getReducedMargin(I);
    I = I.addIndices(Iadd);
end

end