function tests = UT06_Curvature()
tests = functiontests(localfunctions);
end

function testSimulateCurvatureSimple(a)
[m, con, obj, opts] = simple_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

sim = SimulateCurvature(m, con, obj.tF, opts);

d2xdT2 = sim.d2xdT2(4);
a.verifyEqual(numel(d2xdT2), m.nx*nT*nT)
d2udT2 = sim.d2udT2(4);
a.verifyEqual(numel(d2udT2), m.nu*nT*nT)
d2ydT2 = sim.d2ydT2(4);
a.verifyEqual(numel(d2ydT2), m.ny*nT*nT)
end

function testSimulateCurvatureSelectSimple(a)
[m, con, unused, opts] = simple_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

obs = observationSelect([2,4]);
sim = SimulateCurvature(m, con, obs, opts);

d2xdT2 = sim.d2xdT2(:,1);
a.verifyEqual(numel(d2xdT2), m.nx*nT*nT)
d2udT2 = sim.d2udT2(:,1);
a.verifyEqual(numel(d2udT2), m.nu*nT*nT)
d2ydT2 = sim.d2ydT2(:,1);
a.verifyEqual(numel(d2ydT2), m.ny*nT*nT)
end
