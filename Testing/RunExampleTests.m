% Quick and dirty script to run all the example scripts. Just make sure this
% runs w/o errors.

function tests = RunExampleTests()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testT01a(a)
T01a_Building_Model_MassAction;
end

function testT01b(a)
T01b_Building_Model_Analytic;
end

function testT02(a)
T02_Modifying_Models;
end

function testT03a(a)
T03a_Loading_Models_Basic;
end

function testT03b(a)
T03b_Loading_Models_Advanced;
end

function testT04(a)
T04_Basic_ODE_Simulation;
end

function testT05(a)
T05_Fitting_Parameters;
end

function testT06(a)
% empty file
T06_Parameter_Uncertainty;
end

function testT07(a)
% broken, hasn't been updated
T07_Topology_Uncertainty;
end