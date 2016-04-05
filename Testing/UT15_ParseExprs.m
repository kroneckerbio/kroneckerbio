function tests = UT15_ParseExprs()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
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

function testName2ID(a)
if ~verLessThan('matlab', '9.0'); st = warning('off', 'symbolic:sym:sym:DeprecateExpressions'); end

[names, ids, xuvNames, values] = sampleData;
sids = sym(ids);

% Test replacements here
a1 = 'A + B';
b1 = name2id(a1, names, ids, xuvNames);
c1 = values(1) + values(2);
d1 = double(subs(sym(b1), sids, values));
a.verifyEqual(c1 ,d1, 'RelTol', 0.001);

a2 = 'k1*A + k2*B';
b2 = name2id(a2, names, ids, xuvNames);
c2 = values(8)*values(1) + values(9)*values(2);
d2 = double(subs(sym(b2), sids, values));
a.verifyEqual(c2 ,d2, 'RelTol', 0.001);

a3 = 'A_0*B_0 + 2.5*"C:C_0"';
b3 = name2id(a3, names, ids, xuvNames);
c3 = values(18)*values(19) + 2.5*values(21);
d3 = double(subs(sym(b3), sids, values));
a.verifyEqual(c3 ,d3, 'RelTol', 0.001);

a4 = 'k3*A + PI*1.15*exp("C:C_0")'; % Make sure to give pi as PI
b4 = name2id(a4, names, ids, xuvNames);
c4 = values(10)*values(1) + pi*1.15*exp(values(21));
d4 = double(subs(sym(b4), sids, values));
a.verifyEqual(c4 ,d4, 'RelTol', 0.001);

if ~verLessThan('matlab', '9.0') && strcmp(st.state, 'on'); warning('on', 'symbolic:sym:sym:DeprecateExpressions'); end
end

function [names, ids, xuvNames, values] = sampleData
names = {'A' 'B' 'C' 'C:C' 'Bx' 'D' 'cell' 'k1' 'k2' 'k3' 'k4' 'k5' ... % (1:12)
    'k6' 'k7' 'k8' 'k9' 'k10' 'A_0' 'B_0' 'C_0' 'C:C_0' 'Bx_0' 'D_0'}'; % (13:23)
ids = {'mwc479c396_5843_46b6_ad4d_975c12b41199', ...
    'mwa77c4790_692d_4dbe_b4cc_26656a29fd6c', ...
    'mw7002c18e_2352_4b9f_8455_1680c1877b2e', ...
    'mwd9a6e0fe_25c1_4f80_b646_567e7fb71235', ...
    'mwfb297d44_1f12_48c2_a354_8b877faf7847', ...
    'mwe1e554a8_edb7_4f70_a579_47d90c9209c0', ...
    'mw3b5cf017_4336_4325_9585_1a156354bd36', ...
    'mw10cfde8e_4d05_444d_8f86_fb8111c22789', ...
    'mw8ee0260c_16db_41cf_8f15_aba87b282306', ...
    'mw0d002d0c_b4c8_41eb_8c61_e566b9ab042e', ...
    'mw55a40e97_ae47_47ed_a2d8_ebe0e58552df', ...
    'mw4c299a9a_0662_4acf_ad8e_113c22db5e66', ...
    'mwbf857f05_da9c_41e8_90be_aed5cae6a38f', ...
    'mwaa35c74c_3076_4e34_8b55_e766a1f3417c', ...
    'mw8b38bec7_711e_402c_9713_6fb52959ad64', ...
    'mwde2eb888_54ed_40ec_bb95_3dc8d8a5ef5c', ...
    'mw39cb2c47_c9c6_4823_b676_535106b68558', ...
    'bte2b717dc_c361_4da7_a899_8d778e3d4bce', ...
    'btdc491f5c_8152_4c34_a592_35e1ab87f579', ...
    'bt0947af29_7c5b_4f11_ad0f_d4f66104b732', ...
    'bt43a74300_e016_49e8_a1c3_9a3c2c7d2519', ...
    'btb3abaa97_5bf1_4a11_814d_15dcd79d9c5e', ...
    'bte33786d2_c834_45de_b676_fee837fb96ac'}';
xuvNames = {'cell' 'cell' 'cell' 'cell' 'cell' 'cell'}';
values = rand(size(names));

% Sanity checks on sample data
assert(all(size(names)) == all(size(ids)), 'size of names and ids don''t match')
assert(length(xuvNames) < length(names), 'more species volumes than possible species')
assert(all(size(names)) == all(size(values)), 'size of names and values don''t match')
end

function testSimBioExprParser(a)
exprsIn = {'cpt.spc', 'cpt.[spc*]', '[cpt*].spc', '[cpt*].[spc*]', ...
    'cpt.[spc*]+cptcpt.[spc*]', 'cptcpt.[spc*]+cpt.[spc*]'};
exprsExpected = {'"cpt.spc"', '"cpt.spc*"', '"cpt*.spc"', '"cpt*.spc*"', ...
    '"cpt.spc*"+"cptcpt.spc*"', '"cptcpt.spc*"+"cpt.spc*"'};
assert(length(exprsIn) == length(exprsExpected))

for i = 1:length(exprsIn)
    exprOut = simbioExpr2kroneckerbioExpr(exprsIn{i});
    a.verifyEqual(exprOut, exprsExpected{i});
end
end