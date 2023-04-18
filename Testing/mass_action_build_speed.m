function build_times = mass_action_build_speed()
% Builds a mass action model with various numbers of reactions to determine
% how the build speed scales with model size.

n_reactions = [3 10 30 100 300];
build_times = nan(numel(n_reactions), 1);
for i = 1:numel(n_reactions)
    fprintf('%0g reactions...', n_reactions(i))
    
    tic;
    m = InitializeModelMassActionAmount('build speed');
    for j = 1:n_reactions(i)
        % Add compartment
        comp_j = sprintf('v%d',j);
        m = AddCompartment(m, comp_j, 3, rand);
        
        % Add a state with a seed
        state_j = sprintf('x%d',j);
        seed_j = [state_j '_0'];
        m = AddSeed(m, seed_j, rand);
        m = AddState(m, state_j, comp_j, seed_j);
        
        % Add a parameter
        param_j = sprintf('k%d',j);
        m = AddParameter(m, param_j, rand);
        
        % Add an input
        input_j = sprintf('u%d',j);
        m = AddInput(m, input_j, comp_j, rand);
        
        % Add a reaction
        reaction_j = sprintf('r%d',j);
        % Randomly select zero, one, or two reactants and products from the
        % new states and inputs
        choices = {state_j, input_j, ''};
        reactants = choices(randi(3,[1,2]));
        reactants = reactants(~cellfun(@isempty, reactants));
        products = choices(randi(3,[1,2]));
        products = products(~cellfun(@isempty, products));
        m = AddReaction(m, reaction_j, reactants, products, param_j);
        
        % Add an output
        output_j = sprintf('y%d',j);
        m = AddOutput(m, output_j, {state_j rand; input_j rand});
    end
    m = FinalizeModel(m);
    build_times(i) = toc;
    fprintf('done after %0g seconds.\n', build_times(i))
end

end