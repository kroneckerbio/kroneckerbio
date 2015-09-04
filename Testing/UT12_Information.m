function tests = UT12_Information()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveValueSimpleAnalytic(a)
[m, con, obj, opts] = simple_analytic_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveValueMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end