function finalizeModelMassActionMatrix(m)
    
    xnames = strcat({m.States.Compartment}, '.', {m.States.Name})';
    unames = strcat({m.Inputs.Compartment}, '.', {m.Inputs.Name})';
    knames = {m.Parameters.Name}';
    vnames = {m.Compartments.Name}';

    % Set up information on which reactants are multiplied where in
    % a table tracking each reactants' contribution to the rates.
    % The table has one row for each term in the reaction rate, i.e., 
    % a model with two reactions with rates k1*x1 and k2*x1*x2/V would have
    % six rows, two from reaction 1 and four from reaction 2.
    %
    % The table has the following variables: 
    %   reaction_ind
    %       The index of the reaction
    %   product_ind
    %       The index of the product of terms this row belongs to. For
    %       the zeroth-order table, this is the same as reaction_ind. 
    %       As partial derivatives are taken, multiple product terms are
    %       generated for each reaction (i.e., reaction rate k*x generates
    %       k and x from partial derivatives for x and k, respectively).
    %   term_ind
    %       The index of the term in the vector [x; u; k; v; c].
    %   term_exponent
    %       The exponent of the term in the rate equation.
    
    % Track numbers
    nx = numel(xnames);
    nu = numel(unames);
    nk = numel(knames);
    nv = numel(vnames);
    nxuk = nx + nu + nk;
    nxukv = nx + nu + nk + nv;
    
    % Create a map into the indices of the terms
    terminds = struct;
    xukvnames = [xnames; unames; knames; vnames]
    for i = 1:numel(xukvnames)
        terminds.(xukvnames{i}) = i;
    end
    xinds = terminds; % Since states are first in the list, their map is the same
    
    % Calculate compartment size matrix, such that Vxuk*[x;u;k;1] = v
    Vxuk = zeros(numel(m.Compartments), numel([xnames; unames; knames; 1]));
    for i = 1:numel(m.Compartments)
        vi = m.Compartments(i);
        for j = 1:size(vi.Size, 1)
            [name_ij, coefficient_ij] = m.Compartments(i).Size{j,:};
            if isempty(name_ij)
                % Empty name indicates a constant value
                Vxuk(i, end) = coefficient_ij;
            else
                Vxuk(i, terminds.(name_ij)) = coefficient_ij;
            end
        end
    end
    
    % Create the zeroth-order reaction and stoichiometry tables
    Tr = table;
    max_nrows = nr*4; % Maximum of nr*4 terms, if all reactions are bimolecular
    Tr.product_ind = zeros(max_nrows, 1);
    Tr.term_ind = zeros(max_nrows, 1);
    Tr.term_exponent = zeros(max_nrows, 1);
    ti = 0;
    
    % Create stoichiometric table
    Ts = table;
    max_nrows_s = 4*nr;
    Ts.product_ind = zeros(max_nrows_s, 1);
    Ts.dterm_inds = zeros(max_nrows_s, 0);
    Ts.state_coefficient = zeros(max_nrows_s, 1);
    si = 0;
    
    for ri = 1:nr
        r = m.Reactions(ri);
        
        % Add reactant terms
        for rcti = 1:numel(r.Reactants)
            ti = ti + 1;
            Tr.reaction_ind(ti) = ri;
            Tr.term_ind(ti) = terminds.(r.Reactants{rcti});
            Tr.term_exponent(ti) = 1;
        end
        
        % Add rate constant term
        ti = ti + 1;
        Tr.reaction_ind(ti) = ri;
        Tr.term_ind(ti) = terminds.(r.Parameter{1});
        Tr.term_exponent(ti) = 1;
        
        
        % If bimolecular reaction, add volume term
        is_bimolecular = numel(r.Reactants) > 1;
        if is_bimolecular
            ti = ti + 1;
            % Determine which volume to use
            [~, xinds] = ismember(r.Reactants, xnames);
            comps = m.States(xinds).Compartment;
            dims = [comps.Dimension];
            if strcmp(comps(1).Name, comps(2).Name)
                comp_ind = 1;
            elseif dims(1) == dims(2)
                error(['In reaction %d, reactants are from two different compartments with the same dimension.', ...
                'If reactants are in different compartments, they must have different dimensions.'], ri)
            else
                [~, comp_ind] = max(dims);
            end
            Tr.reaction_ind(ti) = ri;
            Tr.term_ind(ti) = terminds.(comps(comp_ind).Name);
            Tr.term_exponent(ti) = -1;
        end
        
        % If there is a multiplier to the rate constant, multiply the stoichiometric coefficients
        mult = r.Parameter{2};
        if isempty(mult)
            mult = 1;
        end
        
        % Add stoichoimetry
        for rcti = 1:numel(r.Reactants)
            si = si + 1;
            Ts.product_ind(si) = ri;
            Ts.state_ind = xinds.(r.Reactants{rcti});
            Ts.state_coefficient = -1*mult;
        end
        for pdi = 1:numel(r.Products)
            si = si + 1;
            Ts.product_ind(si) = ri;
            Ts.state_ind = xinds.(r.Products{pdi});
            Ts.state_coefficient = 1*mult;
        end
    end
    % For the zeroth-order table, product_ind is the same as reaction_ind
    Tr.product_ind = Tr.reaction_ind;
    Tr = Tr(1:ti, :);
    Ts = Ts(1:si, :);
    
    % Also create a table indicating which reactions are added or subtracted from each 
    % state.
    %
    % The table has the following variables:
    %   product_ind
    %   state_ind
    %       The index of the state.
    %   dterm_inds
    %   state_coefficient
    %       The stoichiometric coefficient of the state.
    %
    % This table can be created for derivatives by joining the zeroth order table
    % to the derivative reaction table on reaction_ind.
    
    function [Trd, Tsd] = derivative(Tr, Ts)
        Trd = cell(height(Tr)*4, 1);
        Tsd = cell(height(Ts)*4, 1);
        ti = 0;
        % For each product...
        for ri = 1:max(Tr.product_ind)
            % Find rows of Tr and Ts that are this product
            Tri = Tr(Tr.product_ind == ri, :);
            Tsi = Ts(Ts.product_ind == ri, :);
            
            % Get terms found in this product that aren't constants
            ri_terms = unique(Tri.term_ind);
            
            % For each term...
            for j = 1:numel(ri_terms)
                % Start with the table for the product
                Trij = Tri;
                Tsij = Tsi;
                
                % Decrease the exponents for the current term
                old_exponents = Trij.term_exponent;
                new_exponents = old_exponents;
                is_ti = Trij.term_ind == ri_terms(j);
                assert(sum(is_ti) == 1)
                new_exponents(is_ti) = new_exponents(is_ti) - 1;
                Trij.term_exponent = new_exponents;
                
                % Eliminate rows with exponents of 0
                eliminate_row = new_exponents == 0;
                Trij(eliminate_row,:) = [];
                old_exponents(eliminate_row) = [];
                new_exponents(eliminate_row) = [];
                
                % Multiply by old exponent
                Tsij.state_coefficient = Tsij.state_coefficient .* old_exponents(1);
                
                if ri_terms(j) <= nx + nu + nk % is not v
                    % Note that the derivative was taken with respect to the term
                    Tsij.dterm_inds = cellfun(@(dterm_inds){[dterm_inds, ri_terms(j)]}, Tsij.dterm_inds);
                else % is v
                    % If is a volume term, multiply in the coefficients of the terms
                    % contributing to the volume
                    v_ind = ri_terms - nx - nu - nk;
                    v_terms = find(Vxuk(v_ind,1:end-1)); % Exclude last term because it is the constant term
                    Tsijv = cell(numel(v_terms),1);
                    for vti = 1:numel(v_terms)
                        Tsijv{vti} = Tsij;
                        Tsijv{vti}.state_coefficient = Tsijv{vti}.state_coefficient .* Vxuk(v_ind, v_terms(vti));
                        Tsijv{vti}.dterm_inds = cellfun(@(dterm_inds){[dterm_inds, v_terms(vti)]}, Tsijv{vti}.dterm_inds);
                    end
                    Tsij = vertcat(Tsijv{:});
                end
                
                % Create a new product term for this derivative and store it
                ti = ti + 1;
                Trij.product_ind = repmat(ti, height(Trij), 1);
                Trd{ti} = Trij;
                
                % Store this product/term subtable
                Tsij.product_ind = repmat(ti, height(Tsij), 1);
                Tsd{ti} = Tsij;
            end
        end
        Trd = vertcat(Trd{:});
        Trd.product_ind_max = repmat(ti, height(Trd), 1);
        Tsd = vertcat(Tsd{:});
    end
    
    dat0.Vxuk = Vxuk;
    dat0.nx = nx;
    dat0.nu = nu;
    dat0.nk = nk;
    
    function is_xuk(inds, type)
        switch type
        case 'x'
            val = inds <= nx;
        case 'u'
            val = (inds > nx) & (inds <= nx + nu);
        case 'k'
            val = (inds > nx + nu) & (inds <= nx + nu + nk);
        end
    end
    
    % Generate functions
    d_choices = {'x','u','k'};
    Trd = Tr;
    Tsd = Ts;
    for d_order = 0:2
        % Get the table for this order's derivatives
        if d_order > 0
            [Trd, Tsd] = derivative(Trd, Tsd)
        end
        ncombos = 3^d_order;
        for i = 1:ncombos
            derivs_i = cell(1,d_order);
            [derivs_i{:}] = ind2sub(repmat(3,1,d_order), i);
            derivs_i = [derivs_i{:}];
            derivs_i = d_choices(derivs_i);
            
            % Filter the rows down to appropriate rows for this choice of partial derivatives
            filter_Tsd = true(height(Tsd),1);
            for j = 1:d_order
                filter_Tsd = filter_Tsd & is_xuk(Tsd.dterm_inds(:,j), derivs_i{j});
            end
            Tsdi = Tsdi(filter_Tsd,:);
            Trdi = Trd(ismember(Trd.product_ind, Tsdi.product_ind), :)
            
            deriv_sizes = cellfun(@(n)dat0.(n), derivs_i);
            if d_order == 0
                fun_name = "f";
            else
                fun_name = "d" + string(d_order) + "f" + strjoin(strcat("d", derivs_i(end:-1:1)), "");
                fun_name = regexprep(fun_name, "(d[xuk]){2}", "($1)2")
            end
            m.(fun_name) = generate_function(dat0, Trdi, Tsdi, deriv_sizes);
        end
    end
    
    function f = generate_function(dat0, Tr, Ts, deriv_sizes)
        dat = dat0;
        dat.product_ind = Tr.product_ind;
        dat.term_ind = Tr.term_ind;
        dat.term_has_exponent = Tr.term_exponent ~= 1;
        dat.term_nonunary_exponent = Tr.term_exponent(dat.term_has_exponent);
        dat.product_ind_max = Tr.product_ind_max(1);
        dat.S_product_ind = Ts.product_ind;
        d_order = numel(deriv_sizes);
        if d_order == 0
            dat.nrows = nx;
            dat.ncols = 1;
            row_inds = dat.state_ind;
            col_inds = zeros(size(dat.row_inds,1),0);
        else
            dat.nrows = nx.*nxuk.^(d_order - 1);
            dat.ncols = nxuk;
            dterm_inds_cell = mat2cell(Ts.dterm_inds, size(Ts.dterm_inds,1), ones(size(Ts.dterm_inds,2),1));
            lin_inds = sub2ind([nx, repmat(nxuk, 1, d_order)], Ts.state_ind, dterm_inds_cell{:});
            [row_inds, col_inds] = ind2sub([dat.nrows, dat.ncols], lin_inds);
        end
        dat.S_subs = [row_inds, col_inds];
        dat.state_has_coefficient = Ts.state_coefficient ~= 1;
        dat.state_nonunary_coefficient = Ts.state_coefficient(dat.state_has_coefficient);
        dat.sparse = d_order > 0;
        m.f = @(t,x,u)fun(t, x, u, k, dat);
    end
end

% dfdxu
[r_i, xu_i] = find(P_r);
for i = numel(r_i):-1:1
    P_rx(i, :) = P_r(r_i(i),:);
    P_rx(i, xu_i) = P_rx(i,xu_i) - 1;
    subs_rx(i,)
end
dat.c = P_r;
dat.c(:, u_start:u_end) = 0;
dat.P = P_r - 1;

end

function val = fun(t, x, u, k, dat)
    % dat fields:
    %   .Vxuk
    %   .product_ind
    %   .term_ind
    %   .term_has_exponent
    %   .term_nonunary_exponent
    %   .product_ind_max
    %   .S_product_ind
    %   .S_subs
    %   .state_has_coefficient
    %   .state_nonunary_coefficient
    %   .nx
    %   .ncols
    %   .sparse
    
    % Calculate compartment sizes
    v = dat.Vxuk*[x; u; k; 1];
    
    terms = log([x; u; k; v]);
    terms = terms(dat.term_ind);
    terms(dat.term_has_exponent) = terms(dat.term_has_exponent).*dat.term_nonunary_exponent;
    
    % Calculate products
    P = exp(accumarray(dat.product_ind, terms, [dat.product_ind_max,1]));
    P = P(dat.S_product_ind);
    P(dat.state_has_coefficient) = P(dat.state_has_coefficient).*dat.state_nonunary_coefficient;
    
    val = accumarray(dat.S_subs, P, [dat.nx, dat.ncols], [], [], dat.sparse);
end