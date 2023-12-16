function out = single_step_mut_alt(task_id, N_runs, alpha_star, phi, R0, R, alt_type,mut_size,id_out)
% Performs multiple trials of single-step mutation from randomly assembled
% communities with some alternate conditions for supplementary figs
% alt_type = 0: Tikhonov model with knock outs
% alt_type = 1: Tikhonov model with knock ins
% alt_type = 2: Dirichlet distributed resource vectors w/ partial knockouts
% alt_type = 3: Dirichlet distributed resource vectors w/ global mutations
% If dirichlet mut_size = gamma determines mutation effect size (1 = KO)
% alpha_star = <S*/R>
% phi = <S*/S>

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end

rng(task_id)
organisms = cell(N_runs,2);
harvests = cell(N_runs,2);
deltas = cell(N_runs,2);
parent_inds = zeros(N_runs,1);

lmbda = sqrt(2)*erfinv(1-2*phi);
alpha = alpha_star/phi;

V_tot_root = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
if V_tot_root < 0
    fprintf('Invalid params: decrease alpha_star or phi.\n');
    return
end
epsilon = V_tot_root.^2*(1 - alpha.*I_func(lmbda));
if epsilon < 0
    fprintf('Invalid params: decrease sigma_R_sq, alpha_star or phi.\n');
    return
end
epsilon = sqrt(epsilon);

for i = 1:N_runs
    seed = round(rand*100000000);
    p_val = R0/R;
    R2_val = 0;
    
    p.N = R;
    p.R2 = R2_val;
    p.sparsity = p_val;
    rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
    p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
    while sum(p.capacity < 0) > 0
        fprintf('Negative resource supply generated, retrying.\n');
        p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
    end
    p.alpha = alpha;
    p.seed = round(1/epsilon+seed);
    
    if alt_type == 0 || alt_type == 1
        p.epsilon = epsilon;
        prop = generateOrganismsTikhonov(p);
    elseif alt_type == 2 || alt_type == 3
        p.epsilon = epsilon/R0;
        prop = generateOrganismsDirichlet(p);
    end
    prop.capacity = p.capacity;
    true_deltas = {};
    
    if alt_type == 0 || alt_type == 1
        [harvest, delta] = findEquilibriumTikhonov(prop);
    elseif alt_type == 2 || alt_type == 3
        [harvest, delta] = findEquilibrium(prop);
    end
    
    organisms{i,1} = prop;
    harvests{i,1} = harvest;
    deltas{i,1} = delta;

    % Check for optimization failure
    if delta(1) > 10 
        out = -3
        return
    end

    if alt_type == 0
        threshold = 1e-3*std(1-harvest(1-harvest > 0));
        survivors = delta > -threshold;

        enzymes = prop.enzymesInSpecies(survivors,:);
        inv_harvest = 1./harvest;
        n_lsq = lsqnonneg(1.0*enzymes', inv_harvest');
        n_lsq = n_lsq/sum(n_lsq);

        abundances = zeros(length(deltas),1);
        abundances(survivors) = n_lsq;

        inv_fitness=enzymes.*repmat(reshape(1-harvest,1,length(harvest)),size(enzymes,1),1);
        
        invasion_mat = inv_fitness.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;
        invasion_mat(sum(enzymes,2)==1,:) = 0;

        if max(max(invasion_mat)) == 0
            fprintf('No beneficial mutations, skipping to next community iteration.\n')
            continue
        end
        
        species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
        enz_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));
        

        prop.P = sum(survivors) + 1;
        prop.params.P = sum(survivors) + 1;
        prop.enzymesInSpecies = prop.enzymesInSpecies(survivors,:);
        prop.enzCount = prop.enzCount(survivors);
        prop.costs = prop.costs(survivors);
        parent_inds(i) = species_ind;

        prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
        prop.enzymesInSpecies(end,enz_ind) = 0;
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind)-1];
        prop.costs = [prop.costs; prop.costs(species_ind)-1];
    elseif alt_type == 1
        threshold = 1e-3*std(1-harvest(1-harvest > 0));
        survivors = delta > -threshold;

        enzymes = prop.enzymesInSpecies(survivors,:);
        inv_harvest = 1./harvest;
        n_lsq = lsqnonneg(1.0*enzymes', inv_harvest');
        n_lsq = n_lsq/sum(n_lsq);

        abundances = zeros(length(deltas),1);
        abundances(survivors) = n_lsq;

        enzymes = prop.enzymesInSpecies(survivors,:);

        inv_fitness=(1-enzymes).*repmat(reshape(harvest-1,1,length(harvest)),size(enzymes,1),1);
        
        invasion_mat = inv_fitness.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;

        if max(max(invasion_mat)) == 0
            fprintf('No beneficial mutations, skipping to next community iteration.\n')
            continue
        end
        
        species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
        enz_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));
        
        prop.P = sum(survivors) + 1;
        prop.params.P = sum(survivors) + 1;
        prop.enzymesInSpecies = prop.enzymesInSpecies(survivors,:);
        prop.enzCount = prop.enzCount(survivors);
        prop.costs = prop.costs(survivors);
        parent_inds(i) = species_ind;

        prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
        prop.enzymesInSpecies(end,enz_ind) = 1;
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind)+1];
        prop.costs = [prop.costs; prop.costs(species_ind)+1];
    elseif alt_type == 2
        threshold = 1e-3*std(1-harvest(1-harvest > 0))/R0;
        survivors = delta > -threshold;

        alphas = prop.enzymesInSpecies(survivors,:);
        fitnesses = prop.budget(survivors);
        alphas = alphas./sum(alphas,2);
        inv_harvest = 1./harvest;
        n_lsq = lsqnonneg(alphas', inv_harvest');
        n_lsq = n_lsq.*exp(-fitnesses);
        n_lsq = n_lsq/sum(n_lsq);

        abundances = zeros(length(deltas),1);
        abundances(survivors) = n_lsq;

        enzymes = prop.enzymesInSpecies(survivors,:);

        invasion_mat = 0*enzymes;
        for j = 1:size(enzymes,1)
            enzymes_partial = enzymes(j,:);
            alphas_new = repmat(enzymes_partial,length(enzymes_partial),1);
            alphas_new = alphas_new - eye(length(enzymes_partial))*mut_size/R0;
            alphas_new = alphas_new./sum(alphas_new,2);
            delta_alphas_ko = alphas_new - enzymes_partial;
            inv_fitness=sum(delta_alphas_ko.*harvest,2);
            inv_fitness(enzymes_partial < mut_size/R0) = 0;
            invasion_mat(j,:) = inv_fitness;
        end
        
        invasion_mat = invasion_mat.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;

        if max(max(invasion_mat)) == 0
            fprintf('No beneficial mutations, skipping to next community iteration.\n')
            continue
        end
        
        species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
        enz_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));
        

        prop.P = sum(survivors) + 1;
        prop.params.P = sum(survivors) + 1;
        prop.enzymesInSpecies = prop.enzymesInSpecies(survivors,:);
        prop.enzCount = prop.enzCount(survivors);
        prop.budget = prop.budget(survivors);

        parent_inds(i) = species_ind;
        prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
        prop.enzymesInSpecies(end,enz_ind) = prop.enzymesInSpecies(end,enz_ind) - mut_size/R0;
        prop.enzymesInSpecies(end,:) = prop.enzymesInSpecies(end,:)/sum(prop.enzymesInSpecies(end,:));
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind)];
        prop.budget = [prop.budget; prop.budget(species_ind)];
    elseif alt_type == 3
        threshold = 1e-3*std(1-harvest(1-harvest > 0))/R0;
        survivors = delta > -threshold;

        alphas = prop.enzymesInSpecies(survivors,:);
        fitnesses = prop.budget(survivors);
        alphas = alphas./sum(alphas,2);
        inv_harvest = 1./harvest;
        n_lsq = lsqnonneg(alphas', inv_harvest');
        n_lsq = n_lsq.*exp(-fitnesses);
        n_lsq = n_lsq/sum(n_lsq);

        abundances = zeros(length(deltas),1);
        abundances(survivors) = n_lsq;

        enzymes = prop.enzymesInSpecies(survivors,:);

        invasion_mat = 0*enzymes;
        delta_alphas_mat = zeros(size(enzymes,1),size(enzymes,2),size(enzymes,2));
        for j = 1:size(enzymes,1)
            for k = 1:size(enzymes,2)
                enzymes_partial = enzymes(j,:);
                enzymes2 = enzymes_partial.*normrnd(1,(1/R0)^(1/2),1,R);
                enzymes2(enzymes2 < 0) = 0;
                enzymes2 = enzymes2/sum(enzymes2);
                delta_alphas = enzymes2 - enzymes_partial;
                delta_alphas_mat(j,k,:) = (mut_size/R0)*delta_alphas/norm(delta_alphas);

                inv_fitness=sum(squeeze(delta_alphas_mat(j,k,:)).*harvest');
                invasion_mat(j,k) = inv_fitness;
            end
        end
        
        invasion_mat = invasion_mat.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;

        if max(max(invasion_mat)) == 0
            fprintf('No beneficial mutations, skipping to next community iteration.\n')
            continue
        end
        
        species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
        enz_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));
        

        prop.P = sum(survivors) + 1;
        prop.params.P = sum(survivors) + 1;
        prop.enzymesInSpecies = prop.enzymesInSpecies(survivors,:);
        prop.enzCount = prop.enzCount(survivors);
        prop.budget = prop.budget(survivors);

        parent_inds(i) = species_ind;
        prop.enzymesInSpecies = [prop.enzymesInSpecies; squeeze(prop.enzymesInSpecies(species_ind,:)) + squeeze(delta_alphas_mat(species_ind,enz_ind,:))'];
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind)];
        prop.budget = [prop.budget; prop.budget(species_ind)];
    else
        fprintf('Invalid alt_type\n');
        return
    end

    if alt_type == 0 || alt_type == 1
        [harvest, delta] = findEquilibriumTikhonov(prop);
    elseif alt_type == 2 || alt_type == 3
        [harvest, delta] = findEquilibrium(prop);
    end
    if delta(1) > 10
        out = -4
        return
    end

    organisms{i,2} = prop;
    harvests{i,2} = harvest;
    deltas{i,2} = delta;
end

if alt_type == 0 || alt_type == 1
    fname = "Results/sstep_astar_" + string(alpha_star) + "_phi_" + string(phi) + "_R0_" + string(R0) + "_R_" + string(R) + "_alttype_" + string(alt_type) + "_" + string(task_id - id_out) + ".mat";
elseif alt_type == 2 || alt_type == 3
    fname = "Results/sstep_astar_" + string(alpha_star) + "_phi_" + string(phi) + "_R0_" + string(R0) + "_R_" + string(R) + "_alttype_" + string(alt_type) + "_gamma_" + string(mut_size) + "_" + string(task_id - id_out) + ".mat";

end

%save(fname,'harvests','deltas','parent_inds','organisms','-v7.3')
save(fname,'harvests','deltas','parent_inds','-v7.3')

out = 0;
end
