% This testing function provides a number of tests of the imaginary step
% finite differences calculator of sensitivities, curvatures, gradients,
% and hessians. It is primarily to be used to debug problems with complex
% numbers in Kronecker functions. If you suspect that problems with complex
% numbers are causing tests to fail (such as use of ' instead of .' causing
% unwanted conjugation of complex numbers), use these tests to help
% determine whether the problem can be isolated to the imaginary step
% finite differences routine.
%
% This test isn't run with the normal suite of tests because it is rather
% slow (it takes about 5 minutes to complete) and is somewhat redundant
% (errors that occur here will almost always also occur in the other unit
% tests).
%
% The tests below will often fail because the imaginary step finite
% difference returns slightly different results than the real step finite
% difference. There probably isn't a problem as long as the difference in
% results is small and the comparison of the imaginary step finite
% difference to the analytic integration passes. The failure in this case
% results from roundoff errors in the real step finite difference results
% originating from 1e-8 being too small of a step size. (The imaginary step
% avoids calculating a small difference between two nearly equal large
% numbers and doesn't have roundoff error problems).

function tests = UT00a_ImaginaryStepFD()
tests = functiontests(localfunctions);
end


%% Sensitivity

function tests = testImaginaryStepFD_MassActionSensitivity(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function tests = testImaginaryStepFD_AnalyticSensitivity(a)
[m, con, ~, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);
tGet = 1:10;

verifySensitivity(a, m, con, tGet, opts)
end

function tests = testImaginaryStepFD_MassActionEventSensitivity(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifySensitivityEvent(a, m, con, obs, opts)
end

function tests = testImaginaryStepFD_AnalyticEventSensitivity(a)
[m, con, ~, opts, eve] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);

obs = observationEvents(10, eve);

verifySensitivityEvent(a, m, con, obs, opts)
end

function verifySensitivity(a, m, con, tGet, opts)
obsSelect = observationSelect(tGet);

simint = SimulateSensitivity(m, con, obsSelect, opts);

opts.ImaginaryStep = false;
simreal_simple = FiniteSimulateSensitivity(m, con, obsSelect, opts);

opts.ImaginaryStep = true;
simimag_simple = FiniteSimulateSensitivity(m, con, obsSelect, opts);
% Complex finite sensitivities aren't implemented yet
%simimag_complex = FiniteSimulateSensitivity(m, con, max(tGet), opts);

%a.verifyEqual(simimag_complex.dydT(tGet,1:m.ny), simreal.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(simimag_simple.dydT, simreal_simple.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(simimag_simple.dydT, simint.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')
end

function verifySensitivityEvent(a, m, con, obs, opts)

simint = SimulateSensitivity(m, con, obs, opts);

opts.ImaginaryStep = false;
simreal = FiniteSimulateSensitivity(m, con, obs, opts);

opts.ImaginaryStep = true;
simimag = FiniteSimulateSensitivity(m, con, obs, opts);

a.verifyEqual(simreal.dyedT, simimag.dyedT, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(simint.dyedT, simimag.dyedT, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')

end

%% Curvature

function tests = testImaginaryStepFD_MassActionCurvature(a)

[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifyCurvature(a, m, con, tGet, opts)

end

function tests = testImaginaryStepFD_AnalyticCurvature(a)

[m, con, ~, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);
tGet = 1:10;

verifyCurvature(a, m, con, tGet, opts)

end

function tests = testImaginaryStepFD_MassActionEventCurvature(a)

[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifyCurvatureEvent(a, m, con, obs, opts)

end

function tests = testImaginaryStepFD_AnalyticEventCurvature(a)

[m, con, ~, opts, eve] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);

obs = observationEvents(10, eve);

verifyCurvatureEvent(a, m, con, obs, opts)

end

function verifyCurvature(a, m, con, tGet, opts)

nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);
obsSelect = observationSelect(tGet);

simint = SimulateCurvature(m, con, obsSelect, opts);

opts.ImaginaryStep = false;
simreal_simple = FiniteSimulateCurvature(m, con, obsSelect, opts);

opts.ImaginaryStep = true;
simimag_simple = FiniteSimulateCurvature(m, con, obsSelect, opts);

% Complex finite curvature is not implemented

a.verifyEqual(size(simimag_simple.d2xdT2), [m.nx*nT*nT, numel(tGet)])
a.verifyEqual(size(simimag_simple.d2udT2), [m.nu*nT*nT, numel(tGet)])
a.verifyEqual(size(simimag_simple.d2ydT2), [m.ny*nT*nT, numel(tGet)])

a.verifyEqual(simreal_simple.d2ydT2, simimag_simple.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(simint.d2ydT2, simimag_simple.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')

end

function verifyCurvatureEvent(a, m, con, obs, opts)

simint = SimulateCurvature(m, con, obs, opts);

opts.ImaginaryStep = false;
simreal = FiniteSimulateCurvature(m, con, obs, opts);

opts.ImaginaryStep = true;
simimag = FiniteSimulateCurvature(m, con, obs, opts);

a.verifyEqual(simreal.d2yedT2, simimag.d2yedT2, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(simint.d2yedT2, simimag.d2yedT2, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')

end

%% Gradient

function tests = testImaginaryStepFD_MassActionGradient(a)
[m, con, obj, opts] = simple_model();

verifyGradient(a, m, con, obj, opts)
end

function tests = testImaginaryStepFD_SteadyStateGradient(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);

verifyGradient(a, m, con, obj, opts)
end

function tests = testImaginaryStepFD_AnalyticGradient(a)
[m, con, obj, opts] = michaelis_menten_model();

verifyGradient(a, m, con, obj, opts)
end

function verifyGradient(a, m, con, obj, opts)
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

% Integration of analytic derivatives
Dinteg = ObjectiveGradient(m, con, obj, opts);

% Real FD
opts.ImaginaryStep = false;
Dreal = FiniteObjectiveGradient(m, con, obj, opts);

% Imaginary FD
opts.ImaginaryStep = true;
Dimag = FiniteObjectiveGradient(m, con, obj, opts);

a.verifyEqual(size(Dimag), [nT,1])

a.verifyEqual(Dimag, Dreal, 'RelTol', 0.001, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(Dimag, Dinteg, 'RelTol', 0.001, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')
end

%% Hessian

function tests = testImaginaryStepFD_MassActionHessian(a)
[m, con, obj, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function tests = testImaginaryStepFD_SteadyStateHessian(a)
simpleopts.steadyState = true;
[m, con, obj, opts] = simple_model(simpleopts);
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function tests = testImaginaryStepFD_AnalyticHessian(a)
[m, con, obj, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

verifyHessian(a, m, con, obj, opts)
end

function verifyHessian(a, m, con, obj, opts)
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

% Integration of analytic derivatives
Hinteg = ObjectiveHessian(m, con, obj, opts);

% Finite differences with real step
opts.ImaginaryStep = false;
Hreal = FiniteObjectiveHessian(m, con, obj, opts);

% Finite differences with imaginary step
opts.ImaginaryStep = true;
Himag = FiniteObjectiveHessian(m, con, obj, opts);

a.verifyEqual(size(Himag), [nT,nT])

a.verifyEqual(Hreal, Himag, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step returned different values than using a real step.')
a.verifyEqual(Hinteg, Himag, 'RelTol', 0.001, 'AbsTol', 1e-4, 'Finite differences using an imaginary step differed from integration of analytic derivatives.')
end