% The root testing function for KroneckerBio
function results = RunAllTests(blank)
if nargin < 1
    blank = true;
end

if blank
    clc
end

% Use same random number seed
rng(1)

% Find file names of all unit test files (begin with "UT")
[test_directory, ~, ~] = fileparts(mfilename('fullpath'));
files_in_test_directory = dir(test_directory);
filenames_in_test_directory = vec({files_in_test_directory.name});
filenames_that_are_tests = regexp(filenames_in_test_directory, '^UT[^.]*', 'match', 'once');
test_filenames = filenames_that_are_tests(~cellfun(@isempty, filenames_that_are_tests));
test_functions = cellfun(@str2func, test_filenames, 'UniformOutput', false);

% Gather all tests from test files
unflat_tests = cellfun(@(f)f(), test_functions, 'UniformOutput', false);
tests = [unflat_tests{:}];

% Run tests
results = tests.run;
end
