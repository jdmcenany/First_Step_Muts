function org = generateOrganismsSpecialists(params)
% Adapted from Tikhonov & Monasson 2017
% Generate a community of specialist species with unique resource use

    rng(params.seed);
    params.alpha = 1;
    params.P = round(params.alpha*params.N);
    
    org.N = params.N;
    org.P = params.N;
    org.params = params;
    org.enzymesInSpecies = generateEnzymeTable(params);
    org.enzCount = sum(org.enzymesInSpecies,2);
    org.budget = generateBudgets(params);
end

function budgets = generateBudgets(params)
    xMu = params.epsilon*random('Normal',0,1,[params.P,1]);
    budgets = xMu;
end

function enzymeTbl = generateEnzymeTable(params)
    enzymeTbl = eye(params.P);
end
