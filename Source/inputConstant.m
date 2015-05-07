function inp = inputConstant(m, u)
assert(isnumeric(u) && all(u >= 0) && numel(u) == m.nu, 'KroneckerBio:inputConstant:u', 'u must be a nonegative vector of length m.nu')
u = vec(u);

inp = Input(m, @(t,q)repmat(u,[1,numel(t)]), [], zeros(0,1), @(t,q)zeros(numel(u),numel(q)), @(t,q)zeros(numel(u)*numel(q),numel(q)));
