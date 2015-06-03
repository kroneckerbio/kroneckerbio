function tests = UT12_Information()
tests = functiontests(localfunctions);
end

function testObjectiveValueSimple(a)
[m, con, obj, opts] = simple_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end

function testObjectiveValueMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();

F1 = ObjectiveInformation(m, con, obj, opts);
end