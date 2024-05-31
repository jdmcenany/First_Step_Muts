function org = generateOrganismsDirichlet(params,beta_param)
% Generate organisms with Dirichlet-distributed resource use
% By default shape parameter will be R0/R, but can also be called as input

if nargin < 2
    beta_param = params.sparsity;
end

    rng(params.seed);
    params.P = round(params.alpha*params.N);
    
    org.N = params.N;
    org.P = params.P;
    org.params = params;

    org.enzymesInSpecies = generateEnzymeTable(params, beta_param);
    org.enzCount = 0*sum(org.enzymesInSpecies,2);

    org.budget = generateBudgets(params);
end

function budgets = generateBudgets(params)
    xMu = params.epsilon*random('Normal',0,1,[params.P,1]);
    budgets = xMu;
end

function enzymeTbl = generateEnzymeTable(params, beta_param)
    N = params.N;
    P = params.P;
    enzymeTbl = gamrnd(beta_param*ones(P,N),1);
    enzymeTbl = enzymeTbl./sum(enzymeTbl,2);
end
