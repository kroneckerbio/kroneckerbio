classdef MassActionAmountModel < Model
    
    methods (Access = protected)
        growReactions(this) % overrides Model's use of Rate with Parameter in Reactions
        qualifyNames(this)
        calculateDerivatives(this)
    end
    
    methods
        AddCompartment(this, name, dimension, size)
        AddState(this, name, compartment, ic)
        AddOutput(this, name, expression)
        AddReaction(this, name, reactants, products, kForward, kReverse, compartment)
    end
    
    methods
        function this = MassActionAmountModel(name)
            this@Model(name);
        end
    end
end