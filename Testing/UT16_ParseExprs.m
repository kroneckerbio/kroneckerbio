function tests = UT16_ParseExprs()
tests = functiontests(localfunctions);
end

function testRegExpReplaceSimpleName(a)
% Replace cell.E and "cell.E" (3 subs)
exprsIn = {'cell.E', 'cell.S', 'acell.E', 'cell.Ea', ...
    '"cell.E"', '"cell.S"', '"acell.E"', '"cell.Ea"', ...
    '"a:cell.E"', '"cell.E:a"', '"a:cell.E:a"', ...
    'cell.E+cell.E'};
exprsExpected = {'id', 'cell.S', 'acell.E', 'cell.Ea', ...
    'id', '"cell.S"', '"acell.E"', '"cell.Ea"', ...
    '"a:cell.E"', '"cell.E:a"', '"a:cell.E:a"', ...
    'id+id'};

for i = 1:length(exprsIn)
    expr = replaceSymbolRegex(exprsIn{i}, 'cell.E', 'id');
    a.verifyEqual(expr, exprsExpected{i});
end
end

function testRegExpReplaceComplexName(a)
% Replace "cell.E:S" (2 subs)
exprsIn = {'cell.E', 'cell.S', 'acell.E', 'cell.Ea', ...
    '"cell.E:S"', '"cell.S"', '"acell.E:S"', '"cell.E:Sa"', ...
    '"a:cell.E:S"', '"cell.E:S#a"', '"a#cell.E:S#a"', ...
    '"cell.E:S"+"cell.E:S"'}';
exprsExpected = {'cell.E', 'cell.S', 'acell.E', 'cell.Ea', ...
    'id', '"cell.S"', '"acell.E:S"', '"cell.E:Sa"', ...
    '"a:cell.E:S"', '"cell.E:S#a"', '"a#cell.E:S#a"', ...
    'id+id'}';

for i = 1:length(exprsIn)
    expr = replaceSymbolRegex(exprsIn{i}, '"cell.E:S"', 'id');
    a.verifyEqual(expr, exprsExpected{i});
end
end

function testParseMathWithBasicModel(a)

end