function counts = histcounts(A, bins)

counts = histc(A, bins);

counts = [counts(1:end-2), counts(end-1)+counts(end)];
