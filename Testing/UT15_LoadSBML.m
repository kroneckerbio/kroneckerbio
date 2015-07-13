function tests = UT15_LoadSBML()
tests = functiontests(localfunctions);
end

function testBasicLoadSBML(a)
opts = [];
opts.UseNames = true;
m = LoadModelSbmlAnalytic('test.xml', opts);
verifyLoadedSBMLModel(a, m);
end

% TODO: make this more comprehensive
function verifyLoadedSBMLModel(a, m)
a.verifyEqual(m.add.nv, 1);
a.verifyEqual(m.add.nk, 3);
a.verifyEqual(m.add.ns, 1);
a.verifyEqual(m.add.nu, 3);
a.verifyEqual(m.add.nx, 1);
a.verifyEqual(m.add.nr, 3);
a.verifyEqual(m.add.nz, 0);
end

% Testing adding outputs and verifying results to SBML-imported model
% Note: FinalizeModel is slow so tests have been rolled into 1
function testEnzymeModelOutputsx(a)
m = LoadModelSbmlAnalytic('enzyme-catalysis-basic.xml');

m = AddOutput(m, 'complex', '"E:S"');
m = AddOutput(m, 'product', '"S#P"');
m = AddOutput(m, 'modified_product', '1.5*"S#P"');
m = AddOutput(m, 'random_y', '"E:S" + sqrt(S)');

a.verifyEqual(m.add.ny, 4);
a.verifyEqual(m.add.Outputs(1).Name, 'complex');
a.verifyEqual(m.add.Outputs(1).Expression, '"E:S"');

m = FinalizeModel(m);

% Random values
t = 10*rand;
x = rand(4,1);
u = [];
% Expected output values
y1 = x(2);
y2 = x(4);
y3 = 1.5*x(4);
y4 = x(2) + sqrt(x(3));
yExpected = [y1 y2 y3 y4]';

a.verifyEqual(m.y(t,x,u), yExpected);
end

