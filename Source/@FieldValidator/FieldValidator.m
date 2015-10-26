classdef FieldValidator
    % Contains field (arguments to functions consisting of names, expressions,
    %   etc.) validator functions.
    % Contains the old fix* private functions
    % Note/TODO: Consider making this a package
    
    methods (Static)
        name = ModelName(name)
        
        name = CompartmentName(name)
        size = CompartmentSize(size, dimension)
        volume = CompartmentSizeMassAction(volume, dimension)
        dimension = CompartmentDimension(dimension)
        
        name = ParameterName(name)
        value = ParameterValue(value)
        
        name = SpeciesName(name)
        ic = StateInitialCondition(ic)
        ic = StateInitialConditionMassAction(ic)
        
        default = InputDefaultValue(default)
        
        name = OutputName(name)
        expression = OutputExpression(expression)
        expression = OutputExpressionMassAction(expression)
        
        [name1, name2] = ReactionName(name)
        names = ReactionSpecies(names)
        name = ReactionParameter(name)
        
        expression = RuleExpression(expression)
    end
    
end