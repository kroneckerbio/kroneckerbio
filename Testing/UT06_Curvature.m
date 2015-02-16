function tests = UT06_Curvature()
tests = functiontests(localfunctions);
end

function testSimulateCurvatureSimple(a)
[m, con, unused, opts] = simple_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

sim = SimulateCurvature(m, con, opts);

d2ydT2 = sim.d2ydT2(4);
a.verifyEqual(numel(d2ydT2), m.ny*nT*nT)
end

function testSimulateCurvatureSelectSimple(a)
[m, con, unused, opts] = simple_model();
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);

sim = SimulateCurvatureSelect(m, con, [2,4], opts);

d2ydT2 = sim.d2ydT2(:,1);
a.verifyEqual(numel(d2ydT2), m.ny*nT*nT)
end
