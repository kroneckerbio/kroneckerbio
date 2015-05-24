function dos = doseZero(m)
dos = Dose(m, @(t,h)zeros(m.ns,numel(t)), [], [], @(t,h)zeros(m.ns,0), @(t,h)zeros(0,0));
