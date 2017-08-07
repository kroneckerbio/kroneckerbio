function tests = UT02_LoadingModels()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

% Make sure these tests are updated as the example files change
function testEquilibriumLoading(a)
m = LoadModelMassAction('Equilibrium.txt');
a.verifyEqual(m.Name, 'Equilibrium');
a.verifyEqual({m.Compartments.Name}, {'v'});
a.verifyEqual({m.Parameters.Name}, {'kon','koff'});
a.verifyEqual({m.States.Name}, {'a','b','c'});
a.verifyEqual({m.Reactions.Name}, {'On','Off'});
a.verifyEqual({m.Outputs.Name}, {'A','B','C'});
verifyDerivatives(a, m);
end

function testSimpleLoading(a)
m = LoadModelMassAction('Simple.txt');
a.verifyEqual(m.Name, 'Simple')
verifyDerivatives(a, m);
end

function testTestingLoading(a)
m = LoadModelMassAction('../Tutorial/Testing.txt');
a.verifyEqual(m.Name, 'Testing')
a.verifyEqual({m.Compartments.Name}, {'v1','v2','v3','v4'});
a.verifyEqual({m.Parameters.Name}, {'k1','k2'});
a.verifyEqual({m.Seeds.Name}, {'s1','s2'});
a.verifyEqual({m.Inputs.Name}, {'u4','u5'});
a.verifyEqual({m.States.Name}, {'x1','x2','x3','x4','x5','x6','x7','x8'});
a.verifyEqual({m.Reactions.Name}, {'r1','r2','','','r35','r5','r6','',''});
a.verifyEqual({m.Outputs.Name}, {'y1','y2','y3','y4','y5'});
verifyDerivatives(a, m);
end

function testCompartmentChanging(a)
m = LoadModelMassAction('Equilibrium.txt');
a.verifyEqual(m.v(0, m.x0(m.s), m.u), 1);

m.Compartment(1).Size = 1;
m = FinalizeModel(m);
a.verifyEqual(m.v(0, m.x0(m.s), m.u), 1);
end

function testVariableCompartmentLoading(a)
m = LoadModelMassAction('variable_compartments.txt');
verifyDerivatives(a, m);

v = m.v(0, m.x0(m.s), m.u);
a.verifyEqual(v(1), v(2))
a.verifyEqual(v(1), v(3))

v = m.v(0, m.x0(m.s), [2;3]);
a.verifyEqual(v(2), 2)
a.verifyEqual(v(3), 3)
end

function testChenLoading(a)
m = LoadModelMassAction('Chen.txt');
verifyDerivatives(a, m);
end

%% More comprehensive advanced model loading
function testBasicSbmlLoading(a)
m = LoadModelSbmlAnalytic('test.xml');
m = FinalizeModel(m);

a.verifyEqual(m.nv, 1);
a.verifyEqual(m.nk, 3);
a.verifyEqual(m.ns, 0);
a.verifyEqual(m.nu, 3);
a.verifyEqual(m.nx, 1);
a.verifyEqual(m.nr, 3);
a.verifyEqual(m.nz, 0);
verifyDerivatives(a, m);
end

function testBasicSbmlLoadingWithSeeds(a)
opts = [];
opts.ICsAsSeeds = true;
m = LoadModelSbmlAnalytic('test.xml', opts);
m = FinalizeModel(m);

a.verifyEqual(m.nv, 1);
a.verifyEqual(m.nk, 3);
a.verifyEqual(m.ns, 1);
a.verifyEqual(m.nu, 3);
a.verifyEqual(m.nx, 1);
a.verifyEqual(m.nr, 3);
a.verifyEqual(m.nz, 0);
verifyDerivatives(a, m);
end

function testEnzymeSbmlLoading(a)
m = LoadModelSbmlAnalytic('enzyme-catalysis-basic.xml');

m = AddOutput(m, 'complex', '"E:S"');
m = AddOutput(m, 'product', '"S#P"');
m = AddOutput(m, 'modified_product', '1.5*"S#P"');
m = AddOutput(m, 'random_y', '"E:S" + sqrt(S)');

a.verifyEqual(m.ny, 4);
a.verifyEqual(m.Outputs(1).Name, 'complex');
a.verifyEqual(m.Outputs(1).Expression, '"E:S"');

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

verifyDerivatives(a, m);
end

function testSimpleAnalytic(a)
m = simple_analytic_model();
verifyDerivatives(a, m);
end

function testSimbioSbmlLoading(a)
file = 'simple_analytic.xml';
m1 = LoadModelSbmlAnalytic(file);
m1 = FinalizeModel(m1);
verifyDerivatives(a, m1);

m2 = LoadModelSimBioAnalytic(file);
% TODO: hack simbio to get the correct compartments out
m2.Compartments(3).Dimension = 2;
m2 = FinalizeModel(m2);
verifyDerivatives(a, m2);

compare_analytic_models(a, m1, m2)
end

function testHigherOrderDose(a)
m = higher_order_dose_model();
verifyDerivatives(a, m);

a.verifyEqual(m.x0([4;5]), [16;28;4;20;6])
a.verifyEqual(m.dx0ds([4;5]), sparse([8,0;7,0;0,0;5,4;0,0]))
a.verifyEqual(m.dx0dk([4;5]), sparse([0,0,0;0,0,4;4,0,0;0,0,0;3,2,0]))
a.verifyEqual(m.d2x0ds2([4;5]), sparse([1,9,4], [1,1,2], [2,1,1], 10, 2))
a.verifyEqual(m.d2x0dk2([4;5]), sparse([3,10,5], [1,1,2], [2,1,1], 15, 3))
a.verifyEqual(m.d2x0dkds([4;5]), sparse([2],[3],[1],10,3))
a.verifyEqual(m.d2x0dsdk([4;5]), sparse([12],[1],[1],15,2))
end

function testSimpleMassActionSbmlLoading(a)
warning('off', 'symbolic2massaction:repeatedSpeciesNames'); % suppress warning for repeated species C
m = LoadModelSbmlMassAction('simple_massaction.xml');
warning('on', 'symbolic2massaction:repeatedSpeciesNames'); % reenable warning for continued session

m = FinalizeModel(m);

verifyDerivatives(a, m);
end

function testSimpleMassActionAsAnalyticSbmlLoading(a)
opts = [];
opts.EvaluateExternalFunctions = true; % simple_massaction has x^2 terms, and power needs to be evaluated
m = LoadModelSbmlAnalytic('simple_massaction.xml');
m = FinalizeModel(m, opts);

verifyDerivatives(a, m);
end

function testSimulateSimpleAnalytic(a)
m = simple_analytic_model();

verifyDerivatives(a, m);
end

function testSimpleMassActionSimBiologyLoading(a)
load('simple_massaction_simbio_model.mat')

warning('off', 'symbolic2massaction:repeatedSpeciesNames'); % suppress warning for repeated species C
m = LoadModelSimBioMassAction(simbiomodel);
warning('on', 'symbolic2massaction:repeatedSpeciesNames'); % reenable warning for continued session

m = FinalizeModel(m);

verifyDerivatives(a, m);
end

function verifyDerivatives(a, m)
% Make all the values about the same size
x0 = rand(m.nx,1) + 1;
u0 = rand(m.nu,1) + 1;
k0 = rand(m.nk,1) + 1;
s0 = rand(m.ns,1) + 1;
m = m.Update(k0);

verifyClose(a, x0, @(x)m.f(0,x,u0), @(x)m.dfdx(0,x,u0))
verifyClose(a, u0, @(u)m.f(0,x0,u), @(u)m.dfdu(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'f'), @(k)k_wrapper(k, 'dfdk'))

verifyClose(a, x0, @(x)m.dfdx(0,x,u0), @(x)m.d2fdx2(0,x,u0))
verifyClose(a, u0, @(u)m.dfdu(0,x0,u), @(u)m.d2fdu2(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdk'), @(k)k_wrapper(k, 'd2fdk2'))

verifyClose(a, u0, @(u)m.dfdx(0,x0,u), @(u)m.d2fdudx(0,x0,u))
verifyClose(a, x0, @(x)m.dfdu(0,x,u0), @(x)m.d2fdxdu(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdx'), @(k)k_wrapper(k, 'd2fdkdx'))

verifyClose(a, x0, @(x)m.dfdk(0,x,u0), @(x)m.d2fdxdk(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdu'), @(k)k_wrapper(k, 'd2fdkdu'))
verifyClose(a, u0, @(u)m.dfdk(0,x0,u), @(u)m.d2fdudk(0,x0,u))

verifyClose(a, s0, @(s)m.x0(s), @(s)m.dx0ds(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'x0'), @(k)k_wrapper_x0(k, 'dx0dk'))

verifyClose(a, s0, @(s)m.dx0dk(s), @(s)m.d2x0dsdk(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'dx0dk'), @(k)k_wrapper_x0(k, 'd2x0dk2'))

verifyClose(a, s0, @(s)m.dx0ds(s), @(s)m.d2x0ds2(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'dx0ds'), @(k)k_wrapper_x0(k, 'd2x0dkds'))

    function val = k_wrapper(k, member)
        m_temp = m.Update(k);
        func = m_temp.(member);
        val = func(0, x0, u0);
    end

    function val = k_wrapper_x0(k, member)
        m_temp = m.Update(k);
        func = m_temp.(member);
        val = func(s0);
    end
end

function verifyClose(a, x0, f, dfdx)
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), dfdx_analytic, 'RelTol', 0.001)
end

function compare_analytic_models(a, m1, m2)
a.verifyEqual(m1.Name, m2.Name)
a.verifyEqual(m1.Compartments, m2.Compartments)
a.verifyEqual(m1.Parameters, m2.Parameters)
a.verifyEqual(m1.Seeds, m2.Seeds)
a.verifyEqual(m1.Inputs, m2.Inputs)
a.verifyEqual(m1.States, m2.States)
a.verifyEqual(m1.Reactions, m2.Reactions)
a.verifyEqual(m1.Rules, m2.Rules)
a.verifyEqual(m1.Outputs, m2.Outputs)
end

function testMassActionModelCopy(a)
m = LoadModelMassAction('Equilibrium.txt');
ms = [m; m];

ms(1) = ms(1).Update(ones(m.nk,1));
ms(2) = ms(2).Update(zeros(m.nk,1));

a.verifyEqual(ms(1).k, ones(m.nk,1));
a.verifyEqual(ms(2).k, zeros(m.nk,1));
end

function testAnalyticModelCopy(a)
m = LoadModelSbmlAnalytic('test.xml');
m = FinalizeModel(m);
ms = [m; m];

ms(1) = ms(1).Update(ones(m.nk,1));
ms(2) = ms(2).Update(zeros(m.nk,1));

a.verifyEqual(ms(1).k, ones(m.nk,1));
a.verifyEqual(ms(2).k, zeros(m.nk,1));
end

function testSaveModelFromFiles(a)
m_old = LoadModelMassAction('Equilibrium.txt');
SaveModel(m_old, 'temp.txt');
m_new = LoadModelMassAction('temp.txt');
a.verifyTrue(modelsAreEqual(m_old, m_new))

m_old = LoadModelMassAction('Simple.txt');
SaveModel(m_old, 'temp.txt');
m_new = LoadModelMassAction('temp.txt');
a.verifyTrue(modelsAreEqual(m_old, m_new))

m_old = LoadModelMassAction('DoseModel.txt');
SaveModel(m_old, 'temp.txt');
m_new = LoadModelMassAction('temp.txt');
a.verifyTrue(modelsAreEqual(m_old, m_new))

m_old = LoadModelMassAction('MAPK_DKDP.txt');
SaveModel(m_old, 'temp.txt');
m_new = LoadModelMassAction('temp.txt');
a.verifyTrue(modelsAreEqual(m_old, m_new))

delete('temp.txt')
end

function testSaveModelComponents(a)
% Test the corner cases of mass action components
m = InitializeModelMassActionAmount('my_name');
m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 2, {'x1', 1.2; 'x3', 5.6});
m = AddCompartment(m, 'cell volume', 2, {'x1', 1.2});
m = AddParameter(m, 'k1', 2.0);
m = AddParameter(m, 'parameter 2', 0);
m = AddParameter(m, 'k3', 2);
m = AddSeed(m, 's1', 6.7);
m = AddSeed(m, 'seed 2', 1000);
m = AddInput(m, 'u1', 'v1');
m = AddInput(m, 'u2', 'v2', 4.5);
m = AddInput(m, 'input 1', 'cell volume', 4.5);
m = AddState(m, 'x1', 'v2');
m = AddState(m, 'x2', 'v1', 12.5);
m = AddState(m, 'x3', 'v1', 's1');
m = AddState(m, 'x4', 'v1', {'s1', 4.2});
m = AddState(m, 'state 1', 'cell volume', {'s1', 4.2; '', 2; 'seed 2', 1});
m = AddReaction(m, '', 'x1', 'x3', 'k1');
m = AddReaction(m, '', {'x2', 'x3'}, 'x4', 'parameter 2', 'k1');
m = AddReaction(m, '', {'x2', 'x3'}, {'x4', 'x3', 'x4'}, 'parameter 2');
m = AddReaction(m, 'r 4', {}, 'x1', {'parameter 2', 3}, {'k3', 2});
m = AddOutput(m, 'y1', 'x1');
m = AddOutput(m, 'y2', {'x2', 2; '', 10; 'x3', 4});
m_old = FinalizeModel(m);

SaveModel(m_old, 'temp.txt');
m_new = LoadModelMassAction('temp.txt');
a.verifyTrue(modelsAreEqual(m_old, m_new))

delete('temp.txt')
end
