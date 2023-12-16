function out = multi_step_muts(task_id,N_stop, alpha_star, phi, R0, R, dx_typ,id_out)
% Performs multi-step knockout mutations until one lineage accumulates
% N_step mutations
% alpha_star = <S*/R>
% phi = <S*/S>
% Mutations have a normally distributed change in pure fitness with
% standard deviation dx_typ * sigma_inv

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end
rng(task_id)

N_curr = 0;

lmbda = sqrt(2)*erfinv(1-2*phi);
threshold_prefactor = 1e-3;

alpha = alpha_star/phi;
V_tot_root = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
if V_tot_root < 0
    fprintf('Invalid params: decrease alpha_star or phi.\n');
    return
end
epsilon = V_tot_root*sqrt(1 - alpha.*I_func(lmbda));
dx = dx_typ * V_tot_root*sqrt(alpha*I_func(lmbda))/sqrt(R0*(1-R0/R));

seed = round(rand*100000000);
p_val = R0/R;
R2_val = 0;

p.N = R;
p.R2 = R2_val;
p.sparsity = p_val;
rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N)*rFluctuationMode);
p.alpha = alpha;
p.seed = round(1/epsilon+seed);
p.epsilon = epsilon/R0;
prop = generateOrganisms(p);
prop.capacity = p.capacity;

[harvest1, deltas1] = findEquilibrium(prop);

if deltas1(1) > 10
    out = -3
    return
end

threshold = threshold_prefactor*std(1-harvest1(1-harvest1 > 0))/R0;
survivors = deltas1 > -threshold;

alphas = prop.enzymesInSpecies(survivors,:);
fitnesses = prop.budget(survivors);
alphas = alphas./sum(alphas,2);
inv_harvest = 1./harvest1;
n_lsq = lsqnonneg(alphas', inv_harvest');
n_lsq = n_lsq.*exp(-fitnesses);
n_lsq = n_lsq/sum(n_lsq);

abundances = zeros(length(deltas1),1);
abundances(survivors) = n_lsq;

living_inds = find(survivors);
lineages = 1:length(living_inds);
accumulated_mutations = 0*deltas1;

k = 0;

harvest = harvest1;

organisms_init = prop;
all_harvests = cell(2,1);
all_deltas = cell(2,1);
all_abundances = cell(2,1);
species_present = cell(2,1);
species_info = zeros(length(deltas1),3); % generation started, parent species, enz. mutated

all_harvests{1} = harvest1;
all_deltas{1} = deltas1;
species_present{1} = living_inds;
all_abundances{1} = abundances;

while N_curr < N_stop
    k = k + 1;
    if k == 91
        organisms_final = prop;
    end
    
    prop.P = length(lineages);
    prop.params.P = length(lineages);
    prop.enzymesInSpecies = prop.enzymesInSpecies(living_inds, :);
    if k == 1
        enz1 = prop.enzymesInSpecies;
    end
    prop.enzCount = prop.enzCount(living_inds,:);
    prop.budget = prop.budget(living_inds);


    enzymes = prop.enzymesInSpecies;   
    rand_shift = -exprnd(dx/R0,size(enzymes,1),size(enzymes,2));

    alphas_new = enzymes.*(1./(sum(enzymes,2)-1)).*exp(rand_shift);
    alphas_old = enzymes.*(1./sum(enzymes,2));
    delta_alphas_ko = alphas_new - alphas_old;
    inv_fitness=sum(delta_alphas_ko.*harvest,2);
    inv_fitness = inv_fitness - (delta_alphas_ko + 1./sum(enzymes,2)).*harvest;
   
    invasion_mat = inv_fitness.*abundances(living_inds);
    invasion_mat(invasion_mat < 0) = 0;
    invasion_mat(prop.enzymesInSpecies == 0) = 0;
 
    species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
    enzyme_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));

    accumulated_mutations(lineages(species_ind)) = accumulated_mutations(lineages(species_ind)) + 1;
    N_curr = max(accumulated_mutations);
    
    prop.P = prop.P + 1;
    prop.params.P = prop.params.P + 1;
    prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
    prop.enzymesInSpecies(end,enzyme_ind) = 0;
    prop.enzCount = [prop.enzCount; sum(prop.enzymesInSpecies(end,:))];
    prop.budget = [prop.budget; prop.budget(species_ind)+rand_shift(species_ind,enzyme_ind)];
    lineages = [lineages, lineages(species_ind)];
    
    [harvest_next, deltas] = findEquilibrium(prop);
    if deltas1(1) > 10
        out = -3
        break
    end


    threshold = threshold_prefactor*std(harvest(1-harvest>0))/R0;
    survivors = deltas > -threshold;
    harvest = harvest_next;
    
    alphas = prop.enzymesInSpecies(survivors,:);
    fitnesses = prop.budget(survivors);
    alphas = alphas./sum(alphas,2);
    inv_harvest = 1./harvest;
    n_lsq = lsqnonneg(alphas', inv_harvest');
    n_lsq = n_lsq.*exp(-fitnesses);
    n_lsq = n_lsq/sum(n_lsq);
    
    abundances = zeros(length(deltas),1);
    abundances(survivors) = n_lsq;
    
    living_inds = find(survivors);
    lineages = lineages(living_inds);

    all_harvests{k+1} = harvest;
    all_deltas{k+1} = deltas;
    all_abundances{k+1} = abundances;
   
    next_ind = size(species_info,1)+1;
    last_species_present = species_present{k};
    last_species_present(length(last_species_present)+1) = next_ind;
    species_present{k+1} = last_species_present(survivors);
    species_info(next_ind,1) = k;
    species_info(next_ind,2) = last_species_present(species_ind);
    species_info(next_ind,3) = enzyme_ind;
end

prop.P = length(lineages);
prop.params.P = length(lineages);
prop.enzymesInSpecies = prop.enzymesInSpecies(living_inds, :);
prop.enzCount = prop.enzCount(living_inds,:);
prop.budget = prop.budget(living_inds);

fname = "Results/mstep_astar_" + string(alpha_star) + "_phi_" + string(phi) + "_R0_" + string(R0) + "_R_" + string(R) + "_dx_" + string(dx_typ) + "_" + string(task_id - id_out) + ".mat";

save(fname, 'all_harvests','all_deltas','all_abundances','species_present','species_info','organisms_init','organisms_final','-v7.3')
out = 0;
end
