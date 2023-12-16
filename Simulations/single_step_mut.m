function out = single_step_mut(task_id, N_runs, alpha_star, phi, R0, R, delta_x, mut_type, sigma_R_sq,id_out)
% Performs multiple trials of single-step mutation from randomly assembled
% communities
% Mutation is knockout if mut_type = -1, knockin if mut_type = 0, and a
% knockout with returning species if mut_type = 1
% mut_type > 0
% alpha_star = <S*/R>
% phi = <S*/S>
% delta_x measured in units of sigma_inv

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
epsilon = V_tot_root.^2*(1 - alpha.*I_func(lmbda)) - sigma_R_sq*R0*(1-R0/R)*(1-alpha_star).^2;
if epsilon < 0
    fprintf('Invalid params: decrease sigma_R_sq, alpha_star or phi.\n');
    return
end
epsilon = sqrt(epsilon);

delta_x_old = delta_x;
delta_x = delta_x*V_tot_root*sqrt(alpha*I_func(lmbda))/sqrt(R0*(1-R0/R));

for i = 1:N_runs
    seed = round(rand*100000000);
    p_val = R0/R;
    R2_val = sigma_R_sq*R;
    
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
    p.epsilon = epsilon/R0;
    prop = generateOrganisms(p);
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

    if mut_type == -1
        threshold = 1e-3*std(1-harvest(1-harvest > -delta_x/R0))/R0;
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

        alphas_new = enzymes.*(1./(sum(enzymes,2)-1))*exp(delta_x/R0);
        alphas_old = enzymes.*(1./sum(enzymes,2));
        delta_alphas_ko = alphas_new - alphas_old;
        inv_fitness=sum(delta_alphas_ko.*harvest,2);
        inv_fitness = inv_fitness - (delta_alphas_ko + 1./sum(enzymes,2)).*harvest;
        
        invasion_mat = inv_fitness.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;
        invasion_mat(prop.enzymesInSpecies(survivors,:) == 0) = 0;
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
        prop.budget = prop.budget(survivors);
        parent_inds(i) = species_ind;
        prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
        prop.enzymesInSpecies(end,enz_ind) = 0;
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind) - 1];
        prop.budget = [prop.budget; prop.budget(species_ind)+delta_x/R0];
    elseif mut_type == 0
        threshold = 1e-3*std(1-harvest(1-harvest > -delta_x/R0))/R0;
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

        alphas_new = enzymes.*(1./(sum(enzymes,2)+1))*exp(delta_x/R0);
        alphas_old = enzymes.*(1./sum(enzymes,2));
        delta_alphas_ko = alphas_new - alphas_old;
        inv_fitness=sum(delta_alphas_ko.*harvest,2);
        inv_fitness = inv_fitness + (1./(1+sum(enzymes,2))).*harvest;
        
        invasion_mat = inv_fitness.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;
        invasion_mat(prop.enzymesInSpecies(survivors,:) == 1) = 0;
        invasion_mat(sum(enzymes,2)==R,:) = 0;

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
        prop.enzymesInSpecies(end,enz_ind) = 1;
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind) + 1];
        prop.budget = [prop.budget; prop.budget(species_ind)+delta_x/R0];

    elseif mut_type == 1
        threshold = 1e-3*std(1-harvest(1-harvest > -delta_x/R0))/R0;
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

        alphas_new = enzymes.*(1./(sum(enzymes,2)-1))*exp(delta_x/R0);
        alphas_old = enzymes.*(1./sum(enzymes,2));
        delta_alphas_ko = alphas_new - alphas_old;
        inv_fitness=sum(delta_alphas_ko.*harvest,2);
        inv_fitness = inv_fitness - (delta_alphas_ko + 1./sum(enzymes,2)).*harvest;
        
        invasion_mat = inv_fitness.*abundances(survivors);
        invasion_mat(invasion_mat < 0) = 0;
        invasion_mat(prop.enzymesInSpecies(survivors,:) == 0) = 0;
        invasion_mat(sum(enzymes,2)==1,:) = 0;

        if max(max(invasion_mat)) == 0
            fprintf('No beneficial mutations, skipping to next community iteration.\n')
            continue
        end
        
        species_ind = randsample(size(invasion_mat,1),1,true,sum(invasion_mat,2));
        enz_ind = randsample(size(invasion_mat,2),1,true,invasion_mat(species_ind,:));
        
        prop.P = prop.P + 1;
        prop.params.P = prop.params.P + 1;
        living_inds = find(survivors);
        species_ind = living_inds(species_ind);
        parent_inds(i) = species_ind;

        prop.enzymesInSpecies = [prop.enzymesInSpecies; prop.enzymesInSpecies(species_ind,:)];
        prop.enzymesInSpecies(end,enz_ind) = 0;
        prop.enzCount = [prop.enzCount; prop.enzCount(species_ind) - 1];
        prop.budget = [prop.budget; prop.budget(species_ind)+delta_x/R0];
    else
        fprintf('Invalid mut_type\n');
        return
    end

    [harvest, delta] = findEquilibrium(prop);
    if delta(1) > 10
        out = -4
        return
    end

    organisms{i,2} = prop;
    harvests{i,2} = harvest;
    deltas{i,2} = delta;
end

fname = "Results/sstep_astar_" + string(alpha_star) + "_phi_" + string(phi) + "_R0_" + string(R0) + "_R_" + string(R) + "_dx_" + string(delta_x_old) + "_sigmaR2_" + string(sigma_R_sq) + "_muttype_" + string(mut_type) + "_" + string(task_id - id_out) + ".mat";

%save(fname,'harvests','deltas','parent_inds','organisms','-v7.3')
save(fname,'harvests','deltas','parent_inds','-v7.3')

out = 0;
end
