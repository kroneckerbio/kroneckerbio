function tests = UT00_Helpers()
tests = functiontests(localfunctions);
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
