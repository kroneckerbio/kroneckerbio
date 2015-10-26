classdef AnalyticModel < Model
    
    methods (Access = protected)
        qualifyNames(this)
        calculateDerivatives(this)
    end
    
    methods
        AddCompartment(this, name, dimension, size)
        AddState(this, name, compartment, ic)
        AddOutput(this, name, expression)
        AddReaction(this, name, reactants, products, forward, reverse, compartment)
        AddRule(this, name, target, expression, type)
    end
    
    methods
        function this = AnalyticModel(name)
            this@Model(name);
        end
    end
end