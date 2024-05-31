function out = single_step_mut_gaussian_uptakes(task_id, N_runs, alpha, epsilon, beta_param, mut_size,R,id_out)
% Performs multiple trials of single-step global-effect mutations from random
% communities, with resource usage drawn from a Dirichlet distribution with
% shape parameter beta_param
% alpha = S/R
% Std(pure fitnesses) = epsilon/(beta_param*R)

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end

delta_x = 0;
R0 = beta_param*R;

rng(task_id)
organisms = cell(N_runs,2);
harvests = cell(N_runs,2);
deltas = cell(N_runs,2);
parent_inds = zeros(N_runs,1);

for i = 1:N_runs
    seed = round(rand*100000000);
    p_val = R0/R;
    R2_val = 1*R;
    
    p.N = R;
    p.R2 = R2_val;
    p.sparsity = p_val;
    rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
    p.capacity = 100*ones(p.N,1);
    p.alpha = alpha;
    p.seed = round(1/epsilon+seed);
    p.epsilon = epsilon/R0;
    prop = generateOrganismsDirichlet(p,beta_param);
    prop.capacity = p.capacity;
    true_deltas = {};
    
    [harvest, delta] = findEquilibrium(prop);
    
    organisms{i,1} = prop;
    harvests{i,1} = harvest;
    deltas{i,1} = delta;

    % Check for optimization failure
    if delta(1) > 10 
        out = -3
        return
    end

    threshold = 1e-3*std(1-harvest(1-harvest > -delta_x/R0))/R0;
    survivors = delta > -threshold;

    alphas = prop.enzymesInSpecies(survivors,:);
    fitnesses = prop.budget(survivors);
    alphas = alphas./sum(alphas,2);
    inv_harvest = prop.capacity'./(harvest*mean(prop.capacity));
    n_lsq = lsqnonneg(alphas', inv_harvest');
    n_lsq = n_lsq.*exp(-fitnesses);
    n_lsq = n_lsq/sum(n_lsq);

    abundances = zeros(length(deltas),1);
    abundances(survivors) = n_lsq;

    enzymes = prop.enzymesInSpecies(survivors,:);

    all_dalphas = zeros(size(enzymes,1),size(enzymes,2),size(enzymes,2));
    invasion_mat = 0*enzymes;
    for j = 1:size(enzymes,1)
        enzymes_partial = enzymes(j,:);
        alphas_new = repmat(enzymes_partial,length(enzymes_partial),1);


        for k = 1:size(enzymes,2)
           alpha_mut = gamrnd(beta_param*ones(1,size(enzymes,2)),1);
           alpha_mut = alpha_mut/sum(alpha_mut);
           alpha_mut = alpha_mut - mean(alpha_mut);
           alphas_new(k,:) = mean(alphas_new(k,:)) + sqrt(1-mut_size)*(alphas_new(k,:) - mean(alphas_new(k,:)) + sqrt(mut_size)*alpha_mut);
        end
        alphas_new(alphas_new < 0) = 0;
        alphas_new = alphas_new./sum(alphas_new,2);
        delta_alphas_ko = alphas_new - enzymes_partial;
        all_dalphas(j,:,:) = delta_alphas_ko;
        inv_fitness=sum(delta_alphas_ko.*harvest,2);
        invasion_mat(j,:) = inv_fitness;
    end

    invasion_mat = invasion_mat.*(abundances(survivors).*exp(fitnesses-mean(fitnesses)));
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
    prop.enzymesInSpecies(end,:) = prop.enzymesInSpecies(end,:) + squeeze(all_dalphas(species_ind,enz_ind,:))';
    prop.enzCount = [prop.enzCount; prop.enzCount(species_ind)];
    prop.budget = [prop.budget; prop.budget(species_ind)+delta_x/R0];

    [harvest, delta] = findEquilibrium(prop);
    if delta(1) > 10
        out = -4
        return
    end

    organisms{i,2} = prop;
    harvests{i,2} = harvest;
    deltas{i,2} = delta;
end

running_tot = 0;
for i = 1:N_runs
    running_tot = running_tot + length(deltas{i,2}-1)/R;
end
running_tot/N_runs

fname = "Results/sstep_gauss_uptakes_alpha_" + string(alpha) + "_eps_" + string(epsilon) + "_scale_" + string(beta_param) + "_R_" + string(R) + "_" + string(task_id - id_out) + ".mat";

%save(fname,'harvests','deltas','parent_inds','organisms','-v7.3')
save(fname,'harvests','deltas','parent_inds','-v7.3')

out = 0;
end
