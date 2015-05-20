function tests = UT15_LoadSBML()
tests = functiontests(localfunctions);
end

function testBasicLoadSBML(a)
opts = [];
opts.UseNames = true;
opts.EvaluateExternalFunctions = true;

s = sbml2Symbolic('test.xml', opts);
% Test symbolic -> analytic conversion in another test

verifySymbolicModel(a, s);
end

% TODO: make this more comprehensive
function verifySymbolicModel(a, s)
a.verifyEqual(s.nv, 1);
a.verifyEqual(s.nk, 3);
a.verifyEqual(s.ns, 1);
a.verifyEqual(s.S, sparse([1,-1,-1]));
end