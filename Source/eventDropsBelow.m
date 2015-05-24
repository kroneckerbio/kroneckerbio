function eve = eventDropsBelow(m, output, threshold)

eve = Event(-1, @(t,y)y(output)-threshold);
