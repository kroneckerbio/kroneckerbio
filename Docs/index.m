%% KroneckerBio Documentation
% Welcome to the KroneckerBio documentation. This documentation is prepared
% for the current pre-release branch%${COMMIT}. For additional help, please ask the 
% kroneckerbio-users@googlegroups.com mailing list. For bug reports and
% development-related issues, please use the Github issue tracker at
% <https://github.com/kroneckerbio/kroneckerbio/issues>.
%
% KroneckerBio is a systems biology modeling toolbox for Matlab. It
% provides an easy-to-use programming interface for building, simulating,
% and analyzing ODE models of biological systems. The toolbox incorporates
% numerous methods developed in the Tidor lab at MIT. 
%
% For mass action models, simulations can be run using sparse matrix
% operations, a fact that KroneckerBio exploits. Because Matlab has fast
% sparse matrix algorithms, simulating mass action models and running
% analyses that are dependent on simulations in KroneckerBio is very fast.
%
% Included in KroneckerBio is a rich set of methods for quantifying the
% uncertainty in the parameters and topology of a model and doing optimal
% experimental design to predict which experiments would be best for
% reducing the remaining uncertainty.

%% Download
% The latest stable version is %${VERSION} and can be downloaded at
% <https://github.com/kroneckerbio/kroneckerbio/releases>.
%
% The current development version can be downloaded or
% cloned from the Github repository at
% <https://github.com/kroneckerbio/kroneckerbio>.

%% Installation
% KroneckerBio is written entirely in Matlab. No compilation is necessary.
% Once the toolbox is downloaded, run the |InitKronecker.m| Matlab script.
% This script simply adds several KroneckerBio folders to the Matlab path
% so that the functions are available.

%% Tutorials
% The following tutorials are provided to guide the user through common
% workflows.
% 
%${TUTORIAL_TOC}

%% FAQs
% We'll periodically update this with issues that we've seen that have
% common resolutions.
%
% No questions have been frequently asked yet...

%% Developing kroneckerbio
% Kroneckerbio is licensed under the MIT X11 license. Contribute by cloning
% the repo and submitting a pull request! 
%
% TODO: more detailed information on best practices, contributor copyrights, etc. 

%% Function Help
% The API documentation can be found here.
%
% TODO: either put the list of functions here or link to an auto-generated page
