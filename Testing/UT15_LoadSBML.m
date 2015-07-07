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

% Testing AddOutputsToSymbolic for an SBML imported model
% Test hardcodes the expected symbolic expressions of the outputs - may be
% unnecessarily fragile due to the vagaries of the symbolic solver
function testEnzymeModelOutput_complex(a)
s = sbml2Symbolic('enzyme-catalysis-basic.xml');
s = AddOutputsToSymbolic(s, 'complex', '"E:S"');
a.verifyEqual(char(s.yNames{1}), 'complex');
a.verifyEqual(char(s.y), 'mw430a79b8_4a88_4b05_9a1e_d597ef1c9315');
end

function testEnzymeModelOutput_product(a)
s = sbml2Symbolic('enzyme-catalysis-basic.xml');
s = AddOutputsToSymbolic(s, 'product', '"S#P"');
a.verifyEqual(char(s.yNames{1}), 'product');
a.verifyEqual(char(s.y), 'mwa31ff37a_96ba_44eb_8f34_050d071d8d12');
end

function testEnzymeModelOutput_modified_product(a)
s = sbml2Symbolic('enzyme-catalysis-basic.xml');
s = AddOutputsToSymbolic(s, 'modified_product', '1.5*"S#P"');
a.verifyEqual(char(s.yNames{1}), 'modified_product');
a.verifyEqual(char(s.y), '1.5*mwa31ff37a_96ba_44eb_8f34_050d071d8d12');
end

function testEnzymeModelOutput_random_y(a)
s = sbml2Symbolic('enzyme-catalysis-basic.xml');
s = AddOutputsToSymbolic(s, 'random_y', '"E:S" + sqrt(S)');
a.verifyEqual(char(s.yNames{1}), 'random_y');
a.verifyEqual(char(s.y), 'mw430a79b8_4a88_4b05_9a1e_d597ef1c9315 + mw43a730c4_7e47_41c9_af54_9d05902e78ea^(1/2)');
end
