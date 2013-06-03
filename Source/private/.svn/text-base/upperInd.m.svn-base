function ind = upperInd(m)

ind = zeros((m * (1 + m)) / 2,1);
current = 0;
for i = 1:m
    ind(current+1:current+i) = (i-1) * m + 1 : (i-1) * m + i;
    current = current + i;
end
