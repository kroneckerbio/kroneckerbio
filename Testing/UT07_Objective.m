function tests = UT07_Objective()
tests = functiontests(localfunctions);
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

G = ObjectiveValue(m, con, obj, opts);
end
