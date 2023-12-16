function org = generateOrganismsDirichlet(params)
% Generate organisms with Dirichlet-distributed resource use

    rng(params.seed);
    params.P = round(params.alpha*params.N);
    
    org.N = params.N;
    org.P = params.P;
    org.params = params;

    org.enzymesInSpecies = generateEnzymeTable(params);
    org.enzCount = 0*sum(org.enzymesInSpecies,2);

    org.budget = generateBudgets(params);
end

function budgets = generateBudgets(params)
    xMu = params.epsilon*random('Normal',0,1,[params.P,1]);
    budgets = xMu;
end

function enzymeTbl = generateEnzymeTable(params)
    N = params.N;
    P = params.P;
    enzymeTbl = gamrnd(params.sparsity*ones(P,N),1);
    enzymeTbl = enzymeTbl./sum(enzymeTbl,2);
end
