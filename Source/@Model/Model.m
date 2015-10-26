classdef Model < handle & matlab.mixin.Heterogeneous
    % Generic Model exposing common interface required of all models.
    % TODO: make abstract and specify public interface by making stuff private
    
    properties
        Name
        Compartments
        Parameters
        Seeds
        Inputs
        States
        Reactions
        Outputs
        Rules
        
        m % public for testing/debugging for now
        
        Update % function handle
        Ready
        
        nv = 0
        nk = 0
        ns = 0
        nu = 0
        nx = 0
        nr = 0
        ny = 0
        nz = 0
    end
    
    % Declarations of methods in separate files
    methods (Access = protected)
        growCompartments(this)
        growParameters(this)
        growSeeds(this)
        growInputs(this)
        growStates(this)
        growReactions(this)
        growOutputs(this)
        growRules(this)
%         initialize(this)
        
        qualifyNames(this)
        extractSpecies(this) % rename this
        calculateDerivatives(this)
    end
    
    methods
        AddSeed(this, name, value)
        AddInput(this, name, compartment, default)
        AddParameter(this, name, value)
        Finalize(this)
    end
    
    methods
        function this = Model(name)
            if nargin < 1
                name = [];
            end
            
            this.Name = FieldValidator.ModelName(name);
            
            this.growCompartments;
            this.growParameters;
            this.growSeeds;
            this.growInputs;
            this.growStates;
            this.growReactions;
            this.growOutputs;
            this.growRules;
%             this.initialize; % not needed - look at the implications of this
            
            this.Ready = false;
        end
    end
    
end