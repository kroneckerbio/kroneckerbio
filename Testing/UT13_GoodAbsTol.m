function tests = UT13_GoodAbsTol()
tests = functiontests(localfunctions);
end

function testGoodAbsTolSimple(a)
[m, con, obj, opts] = simple_model();
sd = sdLinear(0.1, 1);

opts.AbsTol = GoodAbsTol(m, con, sd, opts);

sim = SimulateSystem(m, con, obj, opts);
sim = SimulateSensitivity(m, con, obj, opts);
sim = ObjectiveValue(m, con, obj, opts);
opts.UseAdjoint = true;
sim = ObjectiveGradient(m, con, obj, opts);
opts.UseAdjoint = false;
sim = ObjectiveGradient(m, con, obj, opts);
end

function testGoodAbsTolMichaelisMenten(a)
[m, con, obj, opts] = michaelis_menten_model();
sd = sdLinear(0.1, 1);

opts.AbsTol = GoodAbsTol(m, con, sd, opts);

sim = SimulateSystem(m, con, obj, opts);
sim = SimulateSensitivity(m, con, obj, opts);
sim = ObjectiveValue(m, con, obj, opts);
opts.UseAdjoint = true;
sim = ObjectiveGradient(m, con, obj, opts);
opts.UseAdjoint = false;
sim = ObjectiveGradient(m, con, obj, opts);
end