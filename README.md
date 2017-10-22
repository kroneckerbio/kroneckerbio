KroneckerBio
============

KroneckerBio is a systems biology and QSP modeling toolbox for Matlab. It provides an easy-to-use programming interface for building, simulating, and analyzing ODE models of biological systems. The toolbox incorporates numerous methods developed in the Tidor lab at MIT.

For mass action models, simulations can be run using sparse matrix operations, a fact that KroneckerBio exploits. Because Matlab has fast sparse matrix algorithms, simulating mass action models and running analyses that are dependent on simulations in KroneckerBio is very fast. If mass action is not sufficiently expressive for your model, models can also be defined using arbitrary analytic expressions.

Included in KroneckerBio is a rich set of methods for quantifying the uncertainty in the parameters and topology of a model and doing optimal experimental design to predict which experiments would be best for reducing the remaining uncertainty.

Installation
------------

KroneckerBio is written entirely in Matlab. No compilation is necessary. Simply download the entire source of the [most recent stable release](https://github.com/kroneckerbio/kroneckerbio/releases). KroneckerBio is now installed.

To use the library, run the `InitKronecker.m` Matlab script at the start of any session, which modifies the Matlab path to make the KroneckerBio functions available. Like `import` in Python or `library` in R, this only persists through the current Matlab session. KroneckerBio can be permanently imported by [saving the Matlab path](https://www.mathworks.com/help/matlab/ref/savepath.html) if desired.

Getting Started
---------------

The Motivating Example section below shows how to run a single simulation. In the `Tutorial` folder, the `T0x` scripts cover the basics of model building and import, simulation, and fitting.

The help section of each KroneckerBio function generally contains good information about how to use it.

In addition to the tutorial, we recommend the [mailing list](https://groups.google.com/forum/#!forum/kroneckerbio-users) for help in using KroneckerBio. Bugs can also be reported here or in the issues tab of the [GitHub repository](https://github.com/kroneckerbio/kroneckerbio).

Motivating Example
------------------

Below is an example showing the building, simulating, and plotting of the distributed-kinase distributed-phosphatase MAPK model, a simple mass action model.

```matlab
%% Build model
m = InitializeModelMassActionAmount('MAPK-DKDP');

% Compartments organize the states of a model. This model uses a dummy
% compartment "v" into which all states will go, which is fine for small
% models.
m = AddCompartment(m, 'v', 3, 1);

% Seed parameters are used by the initial conditions.
m = AddSeed(m, 'S', 2);
m = AddSeed(m, 'P', 1);

% Input species are externally defined and their amounts are not affected 
% by the reactions. Here, the species "E" always has an amount of 1. The
% experimental conditions can override this with an arbitrary function of
% time.
m = AddInput(m, 'E', 'v', 1);

% State species are controlled by the reactions. Each one exists in a
% particular compartment and has a particular initial condition. In mass
% action models, the initial conditions are restricted to a linear
% combination of seed parameters.
m = AddState(m, 'S', 'v', 'S');
m = AddState(m, 'E:S', 'v');
m = AddState(m, 'M', 'v');
m = AddState(m, 'E:M', 'v');
m = AddState(m, 'D', 'v');
m = AddState(m, 'P', 'v', 'P');
m = AddState(m, 'P:D', 'v');
m = AddState(m, 'P:M', 'v');

% Outputs are the observable states of the model. In mass action models,
% they are restricted to being linear combinations of species.
m = AddOutput(m, 'S', {'S', 'E:S'});
m = AddOutput(m, 'M', {'M', 'E:M', 'P:M'});
m = AddOutput(m, 'D', {'D', 'P:D'});

% Kinetic parameters are used by the reactions.
m = AddParameter(m, 'k1on', 0.02);
m = AddParameter(m, 'k1off', 1);
m = AddParameter(m, 'k1cat', 0.01);
m = AddParameter(m, 'k2on', 0.032);
m = AddParameter(m, 'k2off', 1);
m = AddParameter(m, 'k2cat', 15);
m = AddParameter(m, 'k3on', 0.045);
m = AddParameter(m, 'k3off', 1);
m = AddParameter(m, 'k3cat', 0.092);
m = AddParameter(m, 'k4on', 0.01);
m = AddParameter(m, 'k4off', 1);
m = AddParameter(m, 'k4cat', 0.5);

% Reactions are defined by reactants, products, a forward rate constant,
% and reverse rate constant.
m = AddReaction(m, '', {'E', 'S'}, 'E:S', 'k1on', 'k1off');
m = AddReaction(m, '', 'E:S', {'E', 'M'}, 'k1cat');
m = AddReaction(m, '', {'E', 'M'}, 'E:M', 'k2on', 'k2off');
m = AddReaction(m, '', 'E:M', {'E', 'D'}, 'k2cat');
m = AddReaction(m, '', {'P', 'D'}, 'P:D', 'k3on', 'k3off');
m = AddReaction(m, '', 'P:D', {'P', 'M'}, 'k3cat');
m = AddReaction(m, '', {'P', 'M'}, 'P:M', 'k4on', 'k4off');
m = AddReaction(m, '', 'P:M', {'P', 'S'}, 'k4cat');

% This builds the system of ODEs from the definitions above.
m = FinalizeModel(m);

%% Simulation
% The experimental conditions object defines the initial conditions, the
% inputs, and the doses. Here, the default values on the model are used.
con = experimentInitialValue(m);

% The observation scheme object defines defines what information is stored
% during a simulation. Here, we store everything, which is convinient for
% small models.
obs = observationAll(10000);

% A simulation is a model run under specific experimental conditions
% recorded under a specific observation scheme.
sim = SimulateSystem(m, con, obs);

%% Plotting
% Simulation objects from observationAll are easy to plot.
plot(sim.t, sim.y(sim.t))
legend({m.Outputs.Name})
```