function tests = UT10_Fitting()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testSimpleFitting(a)
[m, con, obj, opts] = simple_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testSimpleAnalyticFitting(a)
[m, con, obj, opts] = simple_analytic_model();
opts.MaxIter = 2;

Gold = ObjectiveValue(m, con, obj, opts);

[m, con] = FitObjective(m, con, obj, opts);

Gnew = ObjectiveValue(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
end

function testMichaelisMentenFitting(a)
[m, con, obj, opts] = michaelis_menten_model();

Gold = ObjectiveValue(m, con, obj, opts);

timer = tic;
[m, con, G, D] = FitObjective(m, con, obj, opts);
time = toc(timer);

Gnew = ObjectiveValue(m, con, obj, opts);
Dnew = ObjectiveGradient(m, con, obj, opts);

a.verifyLessThan(Gnew, Gold)
a.verifyLessThan(time, 12)

a.verifyEqual(G, Gnew, 'RelTol', 1e-4)
a.verifyEqual(D, Dnew, 'RelTol', 1e-4)
end

function testSimpleParallelFitting(a)

fitopts.MaxIter = 2;
nExperiments = 3;
nTotalTimePoints = 15;
testfun = generateTestParallel('simple', fitopts, nExperiments, nTotalTimePoints);
testfun(a)

end

function testMichaelisMentenParallelFitting(a)

fitopts.MaxIter = 2;
nExperiments = 3;
nTotalTimePoints = 15;
testfun = generateTestParallel('michaelis_menten', fitopts, nExperiments, nTotalTimePoints);
testfun(a)

end

%% Test generation function for parallel tests

function testfun = generateTestParallel(model, fitopts, nExperiments, nTotalTimePoints)

% Check for parallel toolbox. If missing, just skip this test
noParallelToolbox = isempty(ver('distcomp'));
if noParallelToolbox
    testfun = @(a)true;
    return
end

switch model
    case 'simple'
        [m, con, obj, opts] = simple_model();
    case 'simple_analytic'
        [m, con, obj, opts] = simple_analytic_model();
    case 'michaelis_menten'
        [m, con, obj, opts] = michaelis_menten_model();
end

% Get randomized experiments and their fitting data
tF = obj.tF;
[obj, opts] = randomExperimentFittingData(m, con, opts, tF,...
    nExperiments, nTotalTimePoints);

opts = mergestruct(opts, fitopts);

testfun = @test;

    function test(a)
        % Tests whether parallel and serial FitObjective give the same
        % solution
        [mfit_parallel, confit_parallel, G_parallel, D_parallel] = ...
            FitObjectiveParallel(m, repmat(con,nExperiments,1), obj, opts);

        [mfit, confit, G, D] = FitObjective(m, repmat(con,nExperiments,1), obj, opts);

        testEqual = @(par, ser, message) a.verifyEqual(par, ser,...
            'AbsTol', 1e-12, 'RelTol', 1e-9, message);

        testEqual(G_parallel, G, 'Objective function values were not equal')
        testEqual(D_parallel, D, 'Gradients were not equal')
        testEqual(mfit_parallel.k, mfit.k, 'k values were not equal')
        for i_con = 1:nExperiments
            testEqual(confit_parallel(i_con).s, confit(i_con).s, ['s values were not equal in experiment ' num2str(i_con)])
            testEqual(confit_parallel(i_con).q, confit(i_con).q, ['q values were not equal in experiment ' num2str(i_con)])
            testEqual(confit_parallel(i_con).h, confit(i_con).h, ['h values were not equal in experiment ' num2str(i_con)])
        end
    end

end