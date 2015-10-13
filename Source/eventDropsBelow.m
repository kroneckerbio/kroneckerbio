function eve = eventDropsBelow(m, output, threshold)

% real() call is needed for imaginary finite differences to work with this
% event
eve = Event(-1, @(t,y)real(y(output))-threshold);
