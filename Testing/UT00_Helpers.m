function tests = UT00_Helpers()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testinsertstruct(a)
s = struct('a', {1,2;3,4}, 'b', {5,6;7,8});

b = struct('a', 9);

t = struct('a', {1,2;9,4}, 'b', {5,6;[],8});
a.verifyEqual(insertstruct(s, b, 2), t)
a.verifyEqual(insertstruct(s, b, 2,1), t)

b = struct('d', 4);

t = struct('a', {1,2;3,[]}, 'b', {5,6;7,[]}, 'd', {[],[];[],4});
a.verifyEqual(insertstruct(s, b, 2,2), t)

b = struct('a', 7, 'd', 9);

t = struct('a', {1,7;3,4}, 'b', {5,[];7,8}, 'd', {[],9;[],[]});
a.verifyEqual(insertstruct(s, b, 1,2), t)

b = struct('a', {7;0}, 'd', {9;4});
t = struct('a', {1,7;3,0}, 'b', {5,[];7,[]}, 'd', {[],9;[],4});
a.verifyEqual(insertstruct(s, b, 1:2,2), t)
end

function testnormbndrnd(a)
% Test a single call (old behavior)
mu = 1;
sigma = 1;
lb = 0.5;
ub = 1.5;
n = 4;
r = normbndrnd(mu, sigma, lb, ub, n);
a.verifyEqual(size(r), [n,1])

% Test sigma = 0
sigma = 0;
r = normbndrnd(mu, sigma, lb, ub, n);
a.verifyEqual(size(r), [n,1])
a.verifyEqual(r, repmat(mu,n,1));

% Test vectorized version with even entries
mu = [0,1];
sigma = [1.5,2.5];
lb = -1;
ub = 3;
n = 10;
r = normbndrnd(mu, sigma, lb, ub, n);
a.verifyEqual(size(r), [n,2])

% Test vectorized version with uneven entries
mu = [0, 1];
sigma = [1.5, 2.5];
lb = -1;
ub = 3;
n = [10, 4];
r = normbndrnd(mu, sigma, lb, ub, n);
a.verifyEqual(size(r), [1,2])
a.verifyEqual(size(r{1}), [n(1),1])
a.verifyEqual(size(r{2}), [n(2),1])
end