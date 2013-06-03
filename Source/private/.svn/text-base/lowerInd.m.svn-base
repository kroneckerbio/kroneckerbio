function ind = lowerInd(m)

ind = zeros((m * (1 + m)) / 2,1);
current = 0;
for i = 1:m
    ind(current+1:current+i) = i : m : (m + 1) * (i - 1) + 1;
    current = current + i;
end
