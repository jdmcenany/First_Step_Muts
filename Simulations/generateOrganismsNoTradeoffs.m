function org = generateOrganismsNoTradeoffs(params)
% Adapted from Tikhonov & Monasson 2017
% Generate organisms with binary resource use where budgets are independent
% of number of resources used

    rng(params.seed);
    params.P = round(params.alpha*params.N);
    
    org.N = params.N;
    org.P = params.P;
    org.params = params;
    org.enzymesInSpecies = generateEnzymeTable(params);
    org.enzCount = sum(org.enzymesInSpecies,2);
    org.budget = generateBudgets(params,org.enzCount);
end

function budgets = generateBudgets(params,enzCount)
    xMu = params.epsilon*random('Normal',0,1,[params.P,1]) + log(enzCount/(params.N*params.sparsity));
    budgets = xMu;
end

function enzymeTbl = generateEnzymeTable(params)
    N = params.N;
    P = params.P;
    % Organisms: each pathway with probability 1/2
    enzymeTbl = rand(P,N)<params.sparsity;

    % Must have no duplicates 
    % !>> Can't check for that! Typical N's are too large. <<!
%     organismId = enzymeTbl*(2.^[0:(N-1)])';
%     dupl = length(unique(organismId))<length(organismId);
%     zero = any(organismId==0);
%     if dupl, fprintf('!! Duplicate organisms found !!\n'), end;
%     if zero, fprintf('!! All-zero organism found !!\n'), end;
%     if dupl || zero
%         fprintf('\t Retrying.\n');
%         enzymeTbl = generateEnzymeTable(params);
%     end
    % Check zeros only 
    zero = any(sum(enzymeTbl,2)==0);
    if zero
        fprintf('\t Empty organism generated. Retrying.\n');
        enzymeTbl = generateEnzymeTable(params);
    end
end
