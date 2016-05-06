function tests = UT12_Information()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testObjectiveInformationSimple(a)
[m, con, obj, opts] = simple_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveInformationSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveInformationMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveInformationCovarianceThetaDependence(a)
% Determine whether ObjectiveInformation can correctly calculate the Fisher
% Information Matrix of a normally distributed output where the mean and
% standard error are both set by a parameter in the model.

stdev = 3;

m = InitializeModelAnalytic();
m = AddParameter(m, 'stdev', stdev);
m = addOutputAnalytic(m, 'stdev', 'stdev');
m = finalizeModelAnalytic(m);

con = experimentInitialValue(m, [], [], [], 'Testing ObjectiveInformation');

obs = observationLinearWeightedSumOfSquares(1, 0, sdLinear(0,1), 'Testing ObjectiveInformation');
obj = obs.Objective(stdev + rand); % Random value added should make no difference

opts.UseParams = 1;

for normalize = [false true]
    opts.Normalized = normalize;
    F = ObjectiveInformation(m, con, obj, opts);
    
    % Value calculated from negative expectation of second derivative of
    % log likelihood function of normal distribution with mean stdev and
    % standard deviation stdev
    Ftheoretical = 1./stdev.^2 + 2./stdev.^2; 
    if opts.Normalized
        Ftheoretical = stdev.'*Ftheoretical*stdev;
    end
    a.verifyEqual(F, Ftheoretical, 'AbsTol', 1e-9, 'RelTol', 1e-6)
end

end