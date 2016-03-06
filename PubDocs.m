function PubDocs(opts)
%PUBDOCS Publish kroneckerbio HTML documentation
%
% Inputs:
%   opts [ struct ]
%       Options struct with the following fields:
%       .OutputDir [ string {'Docs/html'} ]
%           Directory in which to put the tree of HTML documentation
%       .TempDir [ string {'Docs/html_'} ]
%           Directory in which to put temporary files for building the docs
%       .PubRootDir [ string {''} ]
%           Web root directory where HTML docs will be published. By
%           default, this is blank/root, indicating the docs will be
%           accessible right after the URL hostname. Assign a non-blank
%           path if hosting in a subdirectory.
%       .Version [ 'stable' | {'current'} ]
%           Whether the docs will refer to the stable release or the
%           current unstable dev branch. 'stable' should only be used for
%           formal numbered releases. 'current' refers to any other commit.
%       .BuildTutorials [ {true} | false ]
%           Whether to build the tutorials. Warning: they take a while so
%           if you want to iterate quickly, run with this set to false.
%
% Notes:
% Currently, this documentation only applies to the current unstable
%   development tree and not the last stable release. At the next version
%   number release, the main documentation page should point to the stable
%   version.
% Switching between multiple versions of the docs isn't fully implemented
%   yet - the directory structure to handle this is implemented (things go
%   in the 'current') subdir now but the main landing page doesn't have
%   links. Do this on the next full release.
% I put together a very simple templating system where strings in the
%   form %${X} are replaced by strings assigned to the variable X in this
%   function. Hopefully the sequence %${X} is unique enough to not occur
%   normally... Note: this isn't the bash short $X format - it just looks
%   for %${X} exactly.
% char(10) is the interpreted newline '\n'. Used to insert multiline
%   strings into templates.
% Since not everything's implemented yet, make sure to remove the TODO's
%   and comments when they are implemented.
% Be careful of case-insensitivity in Windows. For example, if you make a
%   folder or file with a wrong case letter, you'll be able to access it
%   using a Windows-based webserver but not a Linux-based one.
%
% TODO:
%   Matlab version numbers should be included somewhere


%% Clean up input args
if nargin < 1
    opts = [];
end

current_path = fileparts(mfilename('fullpath'));

opts_ = [];
opts_.OutputDir = 'Docs/html';
opts_.TempDir = 'Docs/html_';
opts_.PubRootDir = '';
opts_.Version = 'current';
opts_.BuildTutorials = true;

opts = mergestruct(opts_, opts);

%% Assemble final directory paths
OutputDir = fullfile(current_path, opts.OutputDir, opts.PubRootDir);
TempDir = fullfile(current_path, opts.TempDir);
TutorialDir = fullfile(current_path, 'tutorial');

%% Version Check
% Matlab version - already taken care of by built-in publish function

% Main kroneckerbio version
VERSION_ = importdata('VERSION');
VERSION = VERSION_{1};

% Commit version, if available
[status, cmdout] = system('git rev-parse HEAD');
if status == 0
    commit = cmdout(1:7);
else
    commit = '';
end
COMMIT = [' (commit ' commit ')'];

% All docs except the landing page go in a subfolder named VERSION
BaseOutputDir = OutputDir;
switch opts.Version
    case 'stable'
        SUBDIR = VERSION;
    case 'current'
        SUBDIR = 'current';
    otherwise
        error('KroneckerBio:PubDocs:InvalidVersionSpec', 'opts.Version only allows values stable or current')
end
OutputDir = fullfile(OutputDir, SUBDIR);

%% Add required directories
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir);
end

if ~exist(TempDir, 'dir')
    mkdir(TempDir);
end

%% Add temporary directory to path since things inside it will be executed
requiredDirs = {TempDir, TutorialDir};
pathCell = regexp(path, pathsep, 'split');

for i = 1:length(requiredDirs)
    requiredDir = requiredDirs{i};
    if ispc  % Windows is not case-sensitive
        onPath = any(strcmpi(requiredDir, pathCell));
    else
        onPath = any(strcmp(requiredDir, pathCell));
    end
    if ~onPath
        addpath(requiredDir);
    end
end

%% Build Tutorials
fprintf('Building tutorials...\n')
pub_opts = [];
pub_opts.outputDir = fullfile(OutputDir, 'tutorial');

files_in_tutorial_directory = dir(TutorialDir);
filenames_in_tutorial_directory = vec({files_in_tutorial_directory.name});
filenames_that_are_tutorials = regexp(filenames_in_tutorial_directory, '^T[0-9][0-9][^.]*', 'match', 'once');
tutorial_filenames = filenames_that_are_tutorials(~cellfun(@isempty, filenames_that_are_tutorials));
nTutorials = length(tutorial_filenames);

% Build tutorials - Warning: takes a while
if opts.BuildTutorials
    for i = 1:nTutorials
        f = tutorial_filenames{i};
        publish(f, pub_opts);
    end
end

% Make tutorials table of contents to insert into landing page
TutorialOutputDir = [SUBDIR '/tutorial']; % must be forward slashes for website
TUTORIAL_TOC = '';
for i = 1:nTutorials
    f = tutorial_filenames{i};
    TUTORIAL_TOC = [TUTORIAL_TOC '% * <' TutorialOutputDir '/' f '.html ' f '>' char(10)];
end

%% Build Function Help
fprintf('Building function help...\n')

% Make functions table of contents to insert into landing page

%% Build landing page
fprintf('Building landing page...\n')

s = [];
s.VERSION = VERSION;
s.COMMIT = COMMIT;
s.TUTORIAL_TOC = TUTORIAL_TOC;
sub_template([current_path '/Docs/index.m'], [current_path, '/Docs/html_/index.m'], s);

pub_opts = [];
pub_opts.outputDir = BaseOutputDir;
publish([current_path, '/Docs/html_/index.m'], pub_opts);

%% Clean up extraneous file handles - run this during exceptions too
% TODO: Can make file accesses safer by using RAII
fclose all;
fprintf('done.\n')
end

function sub_template(filein, fileout, vars)
% Substitute variables into templates
%
% Inputs:
%   filein [ string ]
%       Name of file on which to perform template substitution
%   fileout [ string ]
%       Resulting file with subs performed
%   vars [ struct ]
%       Struct containing variables to substitute as fieldnames and
%       substituted content as the values.

% TODO: validate filein and fileout names
% TODO: validate vars values are strings

fields = fieldnames(vars);
nFields = length(fields);

fin = fopen(filein);
fout = fopen(fileout, 'w');
while ~feof(fin)
   s = fgetl(fin);
   for i = 1:nFields
       field = fields{i};
       new = vars.(field);
       s = strrep(s, ['%${' field '}'], new);
   end
   fprintf(fout, '%s\n', s);
end
fclose(fin);
fclose(fout);

end