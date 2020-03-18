function m = finalizeModelMassActionMatrix(m)
% TODO:
%   - Store names of components instead of indices, then convert to indices
%     at the end, making debugging and understanding the code much easier
%   - Work out a derivation of this approach
%   - Use the derivation to sensibly rename the components
%   - Consider reworking the tables so that each row of Tr is a product
%     instead of a term and each row of Ts is also a product

    xnames = strcat({m.States.Compartment}, '.', {m.States.Name})';
    unames = strcat({m.Inputs.Compartment}, '.', {m.Inputs.Name})';
    knames = {m.Parameters.Name}';
    vnames = {m.Compartments.Name}';
    
    k = [m.Parameters.Value].';

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
    nr = numel(m.Reactions);
    nxuk = nx + nu + nk;
    nxukv = nxuk + nv;
    
    % Create a map into the indices of the terms
    terminds = containers.Map([xnames; unames; knames; vnames], 1:nxukv);
    xinds = containers.Map(xnames, 1:nx);
    
    % Also create a map into the compartments for the states and inputs
    comp_xu_map = containers.Map([xnames; unames], [{m.States.Compartment}, {m.Inputs.Compartment}].');
    
    % Calculate compartment size matrix, such that Vxuk*[x;u;k;1] = v
    Vxuk = zeros(numel(m.Compartments), numel([xnames; unames; knames; 1]));
    for i = 1:numel(m.Compartments)
        vi = m.Compartments(i);
        for jj = 1:size(vi.Size, 1)
            [name_ij, coefficient_ij] = m.Compartments(i).Size{jj,:};
            if isempty(name_ij)
                % Empty name indicates a constant value
                Vxuk(i, end) = coefficient_ij;
            else
                Vxuk(i, terminds(name_ij)) = coefficient_ij;
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
            Tr.term_ind(ti) = terminds(r.Reactants{rcti});
            Tr.term_exponent(ti) = 1;
        end
        
        % Add rate constant term
        ti = ti + 1;
        Tr.reaction_ind(ti) = ri;
        Tr.term_ind(ti) = terminds(r.Parameter{1});
        Tr.term_exponent(ti) = 1;
        
        
        % If bimolecular reaction, add volume term
        is_bimolecular = numel(r.Reactants) > 1;
        if is_bimolecular
            ti = ti + 1;
            % Determine which volume to use
            comps = cellfun(@(x){comp_xu_map(x)}, r.Reactants);
            [~,comp_inds] = ismember(comps, {m.Compartments.Name});
            comps = m.Compartments(comp_inds);
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
            Tr.term_ind(ti) = terminds(comps(comp_ind).Name);
            Tr.term_exponent(ti) = -1;
        end
        
        % If there is a multiplier to the rate constant, multiply the stoichiometric coefficients
        mult = r.Parameter{2};
        if isempty(mult)
            mult = 1;
        end
        
        % Add stoichoimetry
        for rcti = 1:numel(r.Reactants)
            tmp = terminds(r.Reactants{rcti});
            if tmp <= nx
                si = si + 1;
                Ts.product_ind(si) = ri+1;
                Ts.state_ind(si) = tmp;
                Ts.state_coefficient(si) = -1*mult;
            end
        end
        for pdi = 1:numel(r.Products)
            tmp = terminds(r.Products{pdi});
            if tmp <= nx
                si = si + 1;
                Ts.product_ind(si) = ri+1;
                Ts.state_ind(si) = tmp;
                Ts.state_coefficient(si) = 1*mult;
            end
        end
    end
    
    % For the zeroth-order table, product_ind is the same as reaction_ind
    Tr.product_ind = Tr.reaction_ind + 1; % Leave index 1 for the constant term
    Tr = Tr(1:ti, :);
    
    % If the same term shows up in a product more than once, simplify the
    % expression by combining them, adding their exponents together
    [gi, Tr_new] = findgroups(Tr(:, {'reaction_ind', 'product_ind', 'term_ind'}));
    Tr_new.term_exponent = splitapply(@sum, Tr.term_exponent, gi);
    Tr = Tr_new;
    Ts = Ts(1:si, :);
    
    % Only keep unique products
    [Tr, Ts] = unique_product(Tr, Ts);

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
        tti = 0;
        % For each product...
        for ppi = 1:max(Tr.product_ind)
            % Find rows of Tr and Ts that are this product
            Tri = Tr(Tr.product_ind == ppi, :);
            Tsi = Ts(Ts.product_ind == ppi, :);
            
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
                if sum(is_ti) ~= 1
                    error
                end
                new_exponents(is_ti) = new_exponents(is_ti) - 1;
                Trij.term_exponent = new_exponents;
                
                % Multiply stoichiometric coefficient by old exponent
                Tsij.state_coefficient = Tsij.state_coefficient .* old_exponents(is_ti);
                
                % Eliminate rows with exponents of 0
                eliminate_row = new_exponents == 0;
                Trij(eliminate_row,:) = [];
                old_exponents(eliminate_row) = [];
%                 new_exponents(eliminate_row) = [];
                 
                if ri_terms(j) <= nx + nu + nk % is not v
                    % Note that the derivative was taken with respect to the term
                    Tsij.dterm_inds = [Tsij.dterm_inds repmat(ri_terms(j), height(Tsij), 1)];
                else % is v
                    % If is a volume term, apply dvdxuk by multiplying in
                    % the coefficients of the terms contributing to the
                    % volume
                    v_ind = ri_terms(j) - nx - nu - nk;
                    v_terms = find(Vxuk(v_ind,1:end-1)); % Exclude last term because it is the constant term
                    Tsijv = cell(numel(v_terms),1);
                    for vti = 1:numel(v_terms)
                        Tsijv{vti} = Tsij;
                        Tsijv{vti}.state_coefficient = Tsijv{vti}.state_coefficient .* Vxuk(v_ind, v_terms(vti));
                        Tsijv{vti}.dterm_inds = [Tsijv{vti}.dterm_inds, repmat(v_terms(vti), height(Tsijv{vti}), 1)];
                    end
                    Tsij = vertcat(Tsijv{:});
                end
                
                % Create a new product term for this derivative and store it
                tti = tti + 1;
                Trij.product_ind = repmat(tti + 1, height(Trij), 1); % Always leave product_ind = 1 for the constant term
                Trd{tti} = Trij;
                
                % Store this product/term subtable
                if isempty(Trij)
                    % If Trij is empty, there are no terms left to multiply. Use
                    % the constant product index of 1
                    Tsij.product_ind = ones(height(Tsij), 1);
                else
                    Tsij.product_ind = repmat(tti + 1, height(Tsij), 1);
                end
                Tsd{tti} = Tsij;
            end
        end
        Trd = vertcat(Trd{:});
        Trd.product_ind_max = repmat(tti, height(Trd), 1);
        Tsd = vertcat(Tsd{:});
        
         % Only keep unique products
        [Trd, Tsd] = unique_product(Trd, Tsd);
    end
    
    dat0.Vxuk = Vxuk;
    dat0.nx = nx;
    dat0.nu = nu;
    dat0.nk = nk;
    
    function val = is_xuk(inds, type)
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
    funs = struct;
    for d_ord = 0:2
        % Get the table for this order's derivatives
        if d_ord > 0
            [Trd, Tsd] = derivative(Trd, Tsd);
        end
        ncombos = 3^d_ord;
        for i = 1:ncombos
            derivs_i = cell(1,d_ord);
            [derivs_i{:}] = ind2sub(repmat(3,1,d_ord), i);
            derivs_i = [derivs_i{:}];
            derivs_i = d_choices(derivs_i);
            
            % Filter the rows down to appropriate rows for this choice of partial derivatives
            filter_Tsd = true(height(Tsd),1);
            for jj = 1:d_ord
                filter_Tsd = filter_Tsd & is_xuk(Tsd.dterm_inds(:,jj), derivs_i{jj});
            end
            Tsdi = Tsd(filter_Tsd,:);
            Trdi = Trd(ismember(Trd.product_ind, Tsdi.product_ind), :);
            
            % Adjust dterm_inds to reflect only a subset of xuk
            for jj = 1:d_ord
                switch derivs_i{jj}
                    case 'x'
                        index_adjust = 0;
                    case 'u'
                        index_adjust = nx;
                    case 'k'
                        index_adjust = nx + nu;
                end
                Tsdi.dterm_inds(:,jj) = Tsdi.dterm_inds(:,jj) - index_adjust;
            end
            
            deriv_sizes = cellfun(@(n)dat0.("n" + n), derivs_i);
            if d_ord == 0
                fun_name = "f";
            elseif d_ord == 1
                fun_name = "dfd" + derivs_i;
            else
                if derivs_i{1} == derivs_i{2}
                    fun_name = "d2fd" + derivs_i{1} + "2";
                else
                    fun_name = "d2f" + strjoin(strcat("d", derivs_i(end:-1:1)), "");
                end
            end
            funs.(fun_name) = generate_function(dat0, Trdi, Tsdi, deriv_sizes);
        end
    end
    
    function f = generate_function(dat0, Tr, Ts, deriv_sizes)
        dat = dat0;
        dat.product_ind = Tr.product_ind;
        dat.term_ind = Tr.term_ind;
        dat.term_exponent = Tr.term_exponent;
        dat.term_has_exponent = Tr.term_exponent ~= 1;
        dat.term_nonunary_exponent = Tr.term_exponent(dat.term_has_exponent);
        dat.product_ind_max = max(Tr.product_ind);
        dat.S_product_ind = Ts.product_ind;
        d_order = numel(deriv_sizes);
        if d_order == 0
            dat.nrows = nx;
            dat.ncols = 1;
            row_inds = Ts.state_ind;
            col_inds = zeros(size(row_inds,1),0);
        else
            dat.nrows = nx.*prod(deriv_sizes(1:end-1));
            if isempty(dat.nrows)
                dat.nrows = nx;
            end
            dat.ncols = deriv_sizes(end);
            dterm_inds_cell = mat2cell(Ts.dterm_inds, size(Ts.dterm_inds,1), ones(size(Ts.dterm_inds,2),1));
            lin_inds = sub2ind([nx, deriv_sizes], Ts.state_ind, dterm_inds_cell{:});
            [row_inds, col_inds] = ind2sub([dat.nrows, dat.ncols], lin_inds);
        end
        dat.S_subs = [row_inds, col_inds];
        dat.state_coefficient = Ts.state_coefficient;
        dat.state_has_coefficient = Ts.state_coefficient ~= 1;
        dat.state_nonunary_coefficient = Ts.state_coefficient(dat.state_has_coefficient);
        dat.sparse = d_order > 0;
        f = @(t, x, u, k)fun(t, x, u, k, dat);
    end

    m.Update = @(k)Update(m, k, funs);
    m = m.Update(k);
end

function val = fun(~, x, u, k, dat)
    % dat fields:
    %   .Vxuk
    %   .product_ind
    %   .term_ind
    %   .term_has_exponent
    %   .term_nonunary_exponent
    %   .term_exponent
    %   .product_ind_max
    %   .S_product_ind
    %   .S_subs
    %   .state_has_coefficient
    %   .state_nonunary_coefficient
    %   .state_coefficient
    %   .nx
    %   .ncols
    %   .sparse
    
    % Calculate compartment sizes
    v = dat.Vxuk*[x; u; k; 1];
    
    % Calculate products
    terms = log([x; u; k; v]);
    P = exp(accumarray(dat.product_ind, terms(dat.term_ind) .* dat.term_exponent, ...
        [dat.product_ind_max,1]));
    
    % Add products to RHS vector/matrix
    val = accumarray(dat.S_subs, P(dat.S_product_ind) .* dat.state_coefficient, ...
        [dat.nrows, dat.ncols], [], [], dat.sparse);
end

function m = Update(m, k, funs)

fnames = fieldnames(funs);
for i = 1:numel(fnames)
    m.(fnames{i}) = @(t,x,u)funs.(fnames{i})(t,x,u,k);
end
m.Update = @(k)Update(m, k, funs);

end

function [Tr, Ts] = unique_product(Tr, Ts)
Tr = sortrows(Tr, {'product_ind', 'term_ind'}); % Sort rows to ensure terms appear in same order
% Create a table of the old product indices, their associated strings,
% and the new product indices
Tr_prods = table;
[unq_product_inds, ~, product_inds_gi] = unique(Tr.product_ind, 'rows', 'stable');
Tr_prods.prod_str = [""; splitapply(@(term_ind, term_exponent)strjoin(string(term_ind) + "^" + string(term_exponent), "*"), ...
    Tr.term_ind, Tr.term_exponent, product_inds_gi)];
Tr_prods.old_product_ind = [1; unq_product_inds];
[~, ~, Tr_prods.new_product_ind] = unique(Tr_prods.prod_str, 'stable');
Tr_prods.keep_product = false(height(Tr_prods),1);
% Only keep the first instance for each unique product
[~, ia] = unique(Tr_prods.new_product_ind); 
Tr_prods.keep_product(ia) = true;
% Look up Tr product inds in the table
[~, rowinds] = ismember(Tr.product_ind, Tr_prods.old_product_ind);
Tr_keep_product = Tr_prods.keep_product(rowinds);
Tr.product_ind = Tr_prods.new_product_ind(rowinds);
Tr = Tr(Tr_keep_product, :);
[~, rowinds] = ismember(Ts.product_ind, Tr_prods.old_product_ind);
Ts.product_ind = Tr_prods.new_product_ind(rowinds);

% Combine rows of Ts in which the same product_ind is added to the same
% state_ind and dterm_inds. To do so, add the coefficients together
[Ts_new, ~, Ts_gi] = unique(Ts(:,{'product_ind', 'dterm_inds', 'state_ind'}), 'stable');
Ts_new.state_coefficient = splitapply(@sum, Ts.state_coefficient, Ts_gi);
Ts = Ts_new;

% Remove rows with zero coefficients
Ts = Ts(Ts.state_coefficient ~= 0, :);
end