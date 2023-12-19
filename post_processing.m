%% FIGURE 2B: Number of survivors vs. species sampled
clear all

epsilon_vals = [0.001 0.4 1.5 5.5];
alpha_vals = 0.2:0.2:4;
astarvals = zeros(20,4,1000);
inds = zeros(20,4);

R0=40;
R=200;
for i = 1:20
    for i2 = 1:4
        fname = "Raw_Results/Assembly_Sims/alpha_epsilon_params/comm_alpha_" + string(alpha_vals(i)) + "_epsilon_" + string(epsilon_vals(i2)) + "_R0_40_R_200_sigmaR2_0.mat";
        load(fname)
        for j = 1:1000
            h = harvests{j};
            d = deltas{j};
            astarvals(i,i2,j) = sum(d > -1e-3*std(h)/R0)/R;
        end
    end
end

alpha_vals_theory = linspace(0.01,4,1000);
epsilon_vals_theory = logspace(-3,1,1000);
epsilon_vals_theory(1) = 0.001;
epsilon_vals_theory(651) = 0.4;
epsilon_vals_theory(794) = 1.5;
epsilon_vals_theory(935) = 5.5;

astar_vals_theory = zeros(1000);
phi_vals = zeros(1000);
lmbda_init = -9.95;
lmbda_init_last = lmbda_init;

all_lmbdas = [];
for i = 1:length(epsilon_vals_theory)
    lmbda_init = lmbda_init_last;
    for j = 1:length(alpha_vals_theory)
        fun = @(L) lmbda_func(L,alpha_vals_theory(j),epsilon_vals_theory(i),R0,R);
        lmbda = fzero(fun, lmbda_init);
        lmbda_init2 = lmbda_init;
        while j > 2 && lmbda < lmbda_init
            lmbda_init2 = lmbda_init2*1.001;
            lmbda = fzero(fun, lmbda_init2);
        end
        all_lmbdas(end+1) = lmbda;
        astar_vals_theory(i,j) = alpha_vals_theory(j)*erfc(lmbda/sqrt(2))/2;
        phi_vals(i,j) = astar_vals_theory(i,j)/alpha_vals_theory(j);
        lmbda_init = lmbda;
        if j == 1
            lmbda_init_last = lmbda;
        end
    end
end

astar_means = mean(astarvals,3);
astar_stds = std(astarvals,0,3);

save('Processed_Results/Fig2B_data','alpha_vals','astar_means','astar_stds','alpha_vals_theory','astar_vals_theory')

%% FIGURE 2C: Distribution of fitness effects
clear all

R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);


alpha_star = 0.9;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
x_vals = linspace(-0.035,0.035,1000);
S_inv_vals1 = normpdf(x_vals,0, sigma_inv);

alpha_star = 0.5;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
S_inv_vals2 = normpdf(x_vals, 0, sigma_inv);

fitnesses_ko1 = zeros(R^2*1000,1);
ind_ko1 = 1;
fitnesses_ko2 = zeros(R^2*1000,1);
ind_ko2 = 1;
fitnesses_ki1 = zeros(R^2*1000,1);
ind_ki1 = 1;
fitnesses_ki2 = zeros(R^2*1000,1);
ind_ki2 = 1;

fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_0.9_phi_0.1_R0_40_R_200_sigmaR2_0.mat";
load(fname)
for i = 1:1000
    h = harvests{i};
    h = length(h)*h/sum(h);
    enzymes = organisms{i}.enzymesInSpecies;
    survivors = deltas{i} > -1e-3*std(h)/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    delta_alphas_ki = enzymes.*(1./(sum(enzymes,2)+1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    budgets = organisms{i}.budget;
    budgets = budgets(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(budgets);
    overall_shift = -mean(log(overall_shift));
    budgets = budgets + overall_shift;

    for j = 1:size(enzymes,1)
        inv_fitness = sum(delta_alphas_ko(j,:).*h);
        inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
        fitnesses_ko1(ind_ko1:ind_ko1+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
        ind_ko1 = ind_ko1+length(inv_fitness);

        inv_fitness = sum(delta_alphas_ki(j,:).*h);
        inv_fitness = inv_fitness + (1/(1+sum(enzymes(j,:))))*h(~enzymes(j,:));
        fitnesses_ki1(ind_ki1:ind_ki1+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
        ind_ki1 = ind_ki1+length(inv_fitness);
    end
end

fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_0.5_phi_0.1_R0_40_R_200_sigmaR2_0.mat";
load(fname)
for i = 1:1000
    h = harvests{i};
    h = length(h)*h/sum(h);
    enzymes = organisms{i}.enzymesInSpecies;
    survivors = deltas{i} > -1e-3*std(h)/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    delta_alphas_ki = enzymes.*(1./(sum(enzymes,2)+1)-1./sum(enzymes,2));

    alphas = enzymes.*(1./sum(enzymes,2));
    budgets = organisms{i}.budget;
    budgets = budgets(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(budgets);
    overall_shift = -mean(log(overall_shift));
    budgets = budgets + overall_shift;
    for j = 1:size(enzymes,1)
        inv_fitness = sum(delta_alphas_ko(j,:).*h);
        inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
        fitnesses_ko2(ind_ko2:ind_ko2+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
        ind_ko2 = ind_ko2+length(inv_fitness);

        inv_fitness = sum(delta_alphas_ki(j,:).*h);
        inv_fitness = inv_fitness + (1/(1+sum(enzymes(j,:))))*h(~enzymes(j,:));
        fitnesses_ki2(ind_ki2:ind_ki2+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
        ind_ki2 = ind_ki2+length(inv_fitness);
    end
end

edges = (-0.035:0.0001:0.035);
edges2 = (-0.0026:0.0004:0.0026);
bins_ko_1 = zeros(length(edges)-1,1);
bins_ko_2 = zeros(length(edges)-1,1);
bins_ki_1 = zeros(length(edges)-1,1);
bins_ki_2 = zeros(length(edges)-1,1);

figure(5)
h1 = histogram(fitnesses_ko1(1:ind_ko1-1),'Normalization','pdf','BinEdges',edges);

figure(6)
h2 = histogram(fitnesses_ki1(1:ind_ki1-1),'Normalization','pdf','BinEdges',edges);
figure(7)
h3 = histogram(fitnesses_ko2(1:ind_ko2-1),'Normalization','pdf','BinEdges',edges);

figure(8)
h4 = histogram(fitnesses_ki2(1:ind_ki2-1),'Normalization','pdf','BinEdges',edges);

save('Processed_Results/Fig2C_data.mat','edges','h1','h2','h3','h4','x_vals','S_inv_vals1','S_inv_vals2')
close all

%% FIGURE 2D: Sigma_inv vs. niche saturation
clear all

R0 = 10;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.001:0.001:0.999;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv1 = sqrt(q/R)/R0;

R0 = 10;
R = 200;
phi = 0.5;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.001:0.001:0.999;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv2 = sqrt(q/R)/R0;

R0 = 100;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.001:0.001:0.999;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv3 = sqrt(q/R)/R0;

R0 = 100;
R = 200;
phi = 0.5;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.001:0.001:0.999;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv4 = sqrt(q/R)/R0;


all_stds = zeros(19,2,2,1000);
all_inds = ones(19,2,2);
astar_vals = 0.05:0.05:0.95;
phi_vals = [0.1 0.5];

R0 = 10;
for i = 1:19
    for i2 = 1:2
        fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_" + string(astar_vals(i)) + "_phi_" + string(phi_vals(i2)) + "_R0_10_R_200_sigmaR2_0.mat";
        load(fname)
        for i3 = 1:1000
            all_fitnesses = zeros(R^2,1);
            ind = 1;
            h = harvests{i3};
            h = length(h)*h/sum(h);
            enzymes = organisms{i3}.enzymesInSpecies;
            survivors = deltas{i3} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{i3}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
        
            for j = 1:size(enzymes,1)
                if sum(enzymes(j,:)) == 1
                    continue
                end
                inv_fitness = sum(delta_alphas_ko(j,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
                all_fitnesses(ind:ind+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
                ind = ind+length(inv_fitness);
            end

            all_stds(i,1,i2,all_inds(i,1,i2)) = std(all_fitnesses(1:ind-1));
            all_inds(i,1,i2) = all_inds(i,1,i2)+1;
        end
    end
end

R0 = 100;
for i = 1:19
    for i2 = 1:2
        fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_" + string(astar_vals(i)) + "_phi_" + string(phi_vals(i2)) + "_R0_100_R_200_sigmaR2_0.mat";
        load(fname)
        for i3 = 1:1000
            all_fitnesses = zeros(R^2,1);
            ind = 1;
            h = harvests{i3};
            h = length(h)*h/sum(h);
            enzymes = organisms{i3}.enzymesInSpecies;
            survivors = deltas{i3} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{i3}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
        
            for j = 1:size(enzymes,1)
                if sum(enzymes(j,:)) == 1
                    continue
                end
                inv_fitness = sum(delta_alphas_ko(j,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
                all_fitnesses(ind:ind+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));
                ind = ind+length(inv_fitness);
            end

            all_stds(i,2,i2,all_inds(i,2,i2)) = std(all_fitnesses(1:ind-1));
            all_inds(i,2,i2) = all_inds(i,2,i2)+1;
        end
    end
end

sigma_invs = mean(all_stds,4);

save('Processed_Results/Fig2D_data','astar_vals','alpha_star','sigma_invs','sigma_inv1','sigma_inv2','sigma_inv3','sigma_inv4')

%% FIGURE 2E: Monoculture-community fitness correlation vs. niche saturation
clear all

astarvals = 0.1:0.1:0.9;
sigmaR2_vals = [0.007358 0.01213 0.02 0.03297 0.05437];
R0 = 40;
R = 200;
fitnesses_mono = zeros(9,5,R^2*1000);
fitnesses_comm = zeros(9,5,R^2*1000);
inds = ones(9,5);


for i = 1:9
    for i2 = 1:5
        fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_sigmaR2_" + string(sigmaR2_vals(i2)) + ".mat";
        load(fname)
        for i3 = 1:1000
            h = harvests{i3};
            h = length(h)*h/sum(h);
            enzymes = organisms{i3}.enzymesInSpecies;
            survivors = deltas{i3} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{i3}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
        
            for j = 1:size(enzymes,1)
                if sum(enzymes(j,:)) == 1
                    continue
                end
                inv_fitness = sum(delta_alphas_ko(j,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
                fitnesses_comm(i,i2,inds(i,i2):inds(i,i2)+length(inv_fitness)-1) = inv_fitness*exp(budgets(j));

                fitnesses_mono(i,i2,inds(i,i2):inds(i,i2)+length(inv_fitness)-1) = (1-organisms{i3}.capacity(find(enzymes(j,:)))/mean(organisms{i3}.capacity))/sum(enzymes(j,:));
                inds(i,i2) = inds(i,i2)+length(inv_fitness);
            end
        end
    end
end

close all
all_corrs = zeros(9,5);
for i = 1:9
    for i2 = 1:5
        corr = corrcoef(fitnesses_comm(i,i2,1:inds(i,i2)-1),fitnesses_mono(i,i2,1:inds(i,i2)-1));
        all_corrs(i,i2) = corr(1,2);
    end
end


phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.001:0.001:0.999;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;

corr_theory = 1./sqrt(1 + R0^2*sigma_inv.^2./((1-alpha_star).^2.*sigmaR2_vals'));

save('Processed_Results/Fig2E_data','astarvals','alpha_star','all_corrs','corr_theory')

%% FIGURE 3A: P(coex) vs. niche saturation
clear all

coex1 = zeros(14,1);
tot1 = zeros(14,1);
std_h1 = zeros(14,1);

coex2 = zeros(14,1);
tot2 = zeros(14,1);
std_h2 = zeros(14,1);

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];

astar_actual1 = zeros(14,1);
astar_actual2 = zeros(14,1);

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.5_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/200;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        std_h1(i) = std_h1(i) + std(h);
        tot1(i) = tot1(i) + 1;
    end
end

for i = 1:14
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/200;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        std_h2(i) = std_h2(i) + std(h);
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
R = 200;
phi = 0.5;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.001:0.001:0.999;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob1 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.001:0.001:0.999;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob2 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));


save('Processed_Results/Fig3A_data.mat','alpha_star_vals','astar_actual1','astar_actual2', 'cprob1','cprob2', 'astarvals','coex1','tot1','coex2','tot2')

%% FIGURE 3B: P(coex) vs. R0
clear all

coex1 = zeros(9,1);
tot1 = zeros(9,1);

coex2 = zeros(9,1);
tot2 = zeros(9,1);

R0_vals = [10 25 40 70 100 130 160 175 190];

for i = 1:9
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.3_phi_0.1_R0_" + string(R0_vals(i)) + "_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0_vals(i)
            continue
        end
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0_vals(i));
        tot1(i) = tot1(i) + 1;
    end
end

for i = 1:9
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_" + string(R0_vals(i)) + "_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:size(deltas,1)
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0_vals(i)
            continue
        end
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0_vals(i));
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 5:1:195;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.3;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0./R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R./(R0.*(1-R0/R));
sigma_inv = sqrt(q/R)./R0;
Csq = (V_tot./(R0.^2.*(1-R0/R).*sigma_inv)).^2;

cprob1 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));


phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.8;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0./R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R./(R0.*(1-R0/R));
sigma_inv = sqrt(q/R)./R0;
Csq = (V_tot./(R0.^2.*(1-R0/R).*sigma_inv)).^2;

cprob2 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

save('Processed_Results/Fig3B_data.mat','R0', 'cprob1','cprob2','R0_vals','coex1','tot1','coex2','tot2')

%% FIGURE 3C: P(coex) vs. Delta X
clear all

coex1 = zeros(11,1);
tot1 = zeros(11,1);

coex2 = zeros(11,1);
tot2 = zeros(11,1);

deltaxvals = -1.5:0.3:1.5;

R0=40;
for i = 1:11
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.3_phi_0.1_R0_40_R_200_dx_" +string(deltaxvals(i)) + "_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot1(i) = tot1(i) + 1;
    end
end

for i = 1:11
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_40_R_200_dx_" +string(deltaxvals(i)) + "_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.3;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob1 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));


fun = @ (C,dx,sigma) (sigma.*(2.*(-power(exp(1),-power(C - dx,2)./(2..*power(sigma,2))) + power(exp(1),-power(dx,2)./(2.*power(sigma,2)))).*sigma + dx.*sqrt(2.*pi).*erf(dx./(sqrt(2).*sigma)) -... 
       dx.*sqrt(2.*pi).*erf((-C + dx)./(sqrt(2).*sigma))))./...
   (2.*(power(sigma,2)./power(exp(1),power(dx,2)./(2.*power(sigma,2))) + dx.*sqrt(pi./2).*sigma.*(1 + erf(dx./(sqrt(2).*sigma)))));

fun2 = @ (gma, R_0, R, Vtot, lma, n, dx, sigma) (-n - lma).*exp(-n.^2./2)./(exp(-lma.^2./2)-lma.*sqrt(pi./2).*erfc(lma./sqrt(2))) .* fun(gma.^2.*Vtot.*(-n-lma)./(R_0.^2.*(1-R_0./R)), dx, sigma);

dx_vals = (-2:0.01:2)*sigma_inv;
pcoex_vals1 = 0*dx_vals;
for i = 1:length(dx_vals)
    pcoex_vals1(i) = integral(@(x) fun2(1,R0,R,V_tot,lmbda,x,dx_vals(i),sigma_inv),-10,-lmbda);
end


phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.8;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob2 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

fun = @ (C,dx,sigma) (sigma.*(2.*(-power(exp(1),-power(C - dx,2)./(2..*power(sigma,2))) + power(exp(1),-power(dx,2)./(2.*power(sigma,2)))).*sigma + dx.*sqrt(2.*pi).*erf(dx./(sqrt(2).*sigma)) -... 
       dx.*sqrt(2.*pi).*erf((-C + dx)./(sqrt(2).*sigma))))./...
   (2.*(power(sigma,2)./power(exp(1),power(dx,2)./(2.*power(sigma,2))) + dx.*sqrt(pi./2).*sigma.*(1 + erf(dx./(sqrt(2).*sigma)))));

fun2 = @ (gma, R_0, R, Vtot, lma, n, dx, sigma) (-n - lma).*exp(-n.^2./2)./(exp(-lma.^2./2)-lma.*sqrt(pi./2).*erfc(lma./sqrt(2))) .* fun(gma.^2.*Vtot.*(-n-lma)./(R_0.^2.*(1-R_0./R)), dx, sigma);

dx_vals = (-2:0.01:2)*sigma_inv;
pcoex_vals2 = 0*dx_vals;
for i = 1:length(dx_vals)
    pcoex_vals2(i) = integral(@(x) fun2(1,R0,R,V_tot,lmbda,x,dx_vals(i),sigma_inv),-10,-lmbda);
end

save('Processed_Results/Fig3C_data.mat','alpha_star_vals', 'pcoex_vals1','pcoex_vals2','coex1','tot1','coex2','tot2')

%% FIGURE 3C INSET: P(Coex) vs. sinv
clear all

R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);


alpha_star_vals = [0.3 0.8];
alpha_vals = alpha_star_vals/phi;

R0 = 40;
V_tot_vals = (1-R0/R)*(2*(1 - alpha_vals*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
epsilon_vals = V_tot_vals.*sqrt(1 - alpha_vals.*I_func(lmbda));
q_vals1 = (R/R0)*(V_tot_vals.^2 - epsilon_vals.^2)/(1-R0/R);
delta_star_vals = sqrt(q_vals1/R)/R0;

S_inv_all = zeros(5000,2);
coex_all = false(5000,2);

fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.3_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    coex_all(k,1) = (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);

    enzymes = organisms{k,1}.enzymesInSpecies;
    survivors = deltas{k,1} > -1e-3*std(h(1-h>0))/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    costs = organisms{k,1}.budget;
    costs = costs(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(costs);
    overall_shift = -mean(log(overall_shift));
    costs = costs + overall_shift;

    inv_fitness = sum(delta_alphas_ko(parent_inds(k),:).*h);
    enz_ind = find(organisms{k,2}.enzymesInSpecies(parent_inds(k),:) ~= organisms{k,2}.enzymesInSpecies(end,:));
    inv_fitness = inv_fitness - (delta_alphas_ko(parent_inds(k),enz_ind) + 1/(sum(enzymes(parent_inds(k),:))))*h(enz_ind);
    inv_fitness = inv_fitness*exp(costs(parent_inds(k)));

    S_inv_all(k,1) = inv_fitness/delta_star_vals(1);
end

fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    coex_all(k,2) = (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);

    enzymes = organisms{k,1}.enzymesInSpecies;
    survivors = deltas{k,1} > -1e-3*std(h(1-h>0))/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    costs = organisms{k,1}.budget;
    costs = costs(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(costs);
    overall_shift = -mean(log(overall_shift));
    costs = costs + overall_shift;

    inv_fitness = sum(delta_alphas_ko(parent_inds(k),:).*h);
    enz_ind = find(organisms{k,2}.enzymesInSpecies(parent_inds(k),:) ~= organisms{k,2}.enzymesInSpecies(end,:));
    inv_fitness = inv_fitness - (delta_alphas_ko(parent_inds(k),enz_ind) + 1/(sum(enzymes(parent_inds(k),:))))*h(enz_ind);
    inv_fitness = inv_fitness*exp(costs(parent_inds(k)));

    S_inv_all(k,2) = inv_fitness/delta_star_vals(2);
end

edges = 0:0.1:1;
cprobs = zeros(length(edges)-1,2,2);
for i = 2:length(edges)
    mask = S_inv_all(:,1) > edges(i-1) & S_inv_all(:,1) < edges(i);
    cprobs(i-1,1,1) = sum(coex_all(mask,1))/sum(mask);
    cprobs(i-1,1,2) = sqrt(sum(coex_all(mask,1)))/sum(mask);

    mask = S_inv_all(:,2) > edges(i-1) & S_inv_all(:,2) < edges(i);
    cprobs(i-1,2,1) = sum(coex_all(mask,2))/sum(mask);
    cprobs(i-1,2,2) = sqrt(sum(coex_all(mask,2)))/sum(mask);
end

R0 = 40;
R = 200;
alpha_star = 0.3;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha = alpha_star/phi;

V_tot = (1-R0/R).*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot.*sqrt(1 - alpha.*I_func(lmbda));
q = (R./R0).*(V_tot.^2 - epsilon.^2)./(1-R0/R);
sigma_inv = sqrt(q_vals1/R)/R0;
S_inv_vals = 0:0.001:2;
C = S_inv_vals*R0^2*(1-R0/R).*sigma_inv(1)./V_tot;

coex_theory1 = (-2 + power(exp(1),power(C + lmbda,2)/2.).*lmbda.*sqrt(2*pi).*erfc((C + lmbda)/sqrt(2)))./...
   (power(exp(1),(C.*(C + 2*lmbda))/2.)*(-2 + power(exp(1),power(lmbda,2)/2.)*lmbda*sqrt(2*pi)*erfc(lmbda/sqrt(2))));

R0 = 40;
R = 200;
alpha_star = 0.8;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha = alpha_star/phi;

V_tot = (1-R0/R).*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot.*sqrt(1 - alpha.*I_func(lmbda));
q = (R./R0).*(V_tot.^2 - epsilon.^2)./(1-R0/R);
S_inv_vals = 0:0.001:2;
C = S_inv_vals*R0^2*(1-R0/R).*sigma_inv(2)./V_tot;

coex_theory2 = (-2 + power(exp(1),power(C + lmbda,2)/2.).*lmbda.*sqrt(2*pi).*erfc((C + lmbda)/sqrt(2)))./...
   (power(exp(1),(C.*(C + 2*lmbda))/2.)*(-2 + power(exp(1),power(lmbda,2)/2.)*lmbda*sqrt(2*pi)*erfc(lmbda/sqrt(2))));

x_vals1 = zeros(5,1);
x_vals2 = zeros(5,1);
x_vals3 = zeros(5,1);
for i = 2:6
    mask = (S_inv_all(:,1) > edges(i-1) & S_inv_all(:,1) < edges(i)) & coex_all(:,1);
    x_vals1(i-1) = mean(S_inv_all(mask,1));

    mask = (S_inv_all(:,2) > edges(i-1) & S_inv_all(:,2) < edges(i)) & coex_all(:,2);
    x_vals2(i-1) = mean(S_inv_all(mask,2));
end

save('Processed_Results/Fig3C_inset_data.mat','S_inv_vals','coex_theory1','coex_theory2','x_vals1','x_vals2','cprobs')

%% FIGURE 4B (with inset): Number of extinctions vs. niche saturation
clear all

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];
num_extinctions_all = zeros(14,10000,2);
num_species_all = zeros(14,10000,2);
tot_all = zeros(14,1);

astar_actual = zeros(14,1);
R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        num_extinctions = sum(d < -1e-3*std(h(1-h>0))/R0) - (d(parent_inds(k)) < -1e-3*std(h(1-h>0))/R0);
        tot_all(i) = tot_all(i)+1;
        num_extinctions_all(i,tot_all(i),1) = num_extinctions;
        num_species_all(i,tot_all(i),1) = length(d)-2;
    end
end

tot_all = zeros(14,1);
R0=16;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_16_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        num_extinctions = sum(d < -1e-3*std(h(1-h>0))/R0) - (d(parent_inds(k)) < -1e-3*std(h(1-h>0))/R0);
        tot_all(i) = tot_all(i)+1;
        num_extinctions_all(i,tot_all(i),2) = num_extinctions;
        num_species_all(i,tot_all(i),2) = length(d)-2;
    end
end

p_extinct = zeros(14,2);
p_extinct_err = zeros(14,2);
for i = 1:14
    for j = 1:2
        p_extinct(i,j) = sum(num_extinctions_all(i,:,j))/sum(num_species_all(i,:,j));
        p_extinct_err(i,j) = sqrt(sum(num_extinctions_all(i,:,j)))/sum(num_species_all(i,:,j));
    end
end


h1=histogram(num_extinctions_all(9,:,1),'BinEdges',0:1:max(num_extinctions_all(9,:,1)+1));

extinction_x = h1.BinEdges(1:end-1);
h1_vals = h1.Values;
poiss_param = poissfit(num_extinctions_all(9,:,1));

save("Processed_Results/Fig4B_data.mat",'astarvals','p_extinct','p_extinct_err','extinction_x','h1_vals','poiss_param','num_species_all')
close all

%% FIGURE 4C: Number of resources shared between mutant and extinct strains
clear all

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];
num_shared_resources = zeros(10000,1);
typical_shared_resources = zeros(10000*500,1);
num_shared_resources2 = zeros(10000,1);
typical_shared_resources2 = zeros(10000*500,1);
tot_all = 0;
tot2 = 1;

R0=40;
fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(9)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    extinctions = find(d < -1e-3*std(h(1-h>0))/R0 & (1:length(d))' ~= parent_inds(k));
    for i = 1:length(extinctions)
        ind = extinctions(i);
        tot_all = tot_all+1;
        num_shared_resources(tot_all) = sum(organisms{k,2}.enzymesInSpecies(end,:) & organisms{k,2}.enzymesInSpecies(ind,:));

        other_enz = organisms{k,2}.enzymesInSpecies(1:end-1,:);
        other_enz(parent_inds(k),:) = [];
        shared_resources = sum(organisms{k,2}.enzymesInSpecies(end,:).*other_enz,2);
        typical_shared_resources(tot2:tot2+length(shared_resources)-1) = shared_resources;
        tot2 = tot2 + length(shared_resources);
    end
end


tot_all = 0;
tot2=1;
R0 = 16;
fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(9)) + "_phi_0.1_R0_16_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    extinctions = find(d < -1e-3*std(h(1-h>0))/R0 & (1:length(d))' ~= parent_inds(k));
    for i = 1:length(extinctions)
        ind = extinctions(i);
        tot_all = tot_all+1;
        num_shared_resources2(tot_all) = sum(organisms{k,2}.enzymesInSpecies(end,:) & organisms{k,2}.enzymesInSpecies(ind,:));
        other_enz = organisms{k,2}.enzymesInSpecies(1:end-1,:);
        other_enz(parent_inds(k),:) = [];
        shared_resources = sum(organisms{k,2}.enzymesInSpecies(end,:).*other_enz,2);            typical_shared_resources2(tot2:tot2+length(shared_resources)-1) = shared_resources;
        tot2 = tot2 + length(shared_resources);
    end
end

edges = 0:(max(num_shared_resources)+1);
h1 = histogram(num_shared_resources,'BinEdges',edges);
figure()
h2 = histogram(num_shared_resources2,'BinEdges',edges);
figure()
h3 = histogram(typical_shared_resources,'BinEdges',edges);
figure()
h4 = histogram(typical_shared_resources2,'BinEdges',edges);
figure()
x_bins = edges(1:end-1);

h1_vals = h1.Values;
h2_vals = h2.Values;
h3_vals = h3.Values;
h4_vals = h4.Values;
overall_mean1 = mean(typical_shared_resources);
overall_mean2 = mean(typical_shared_resources2);

save("Processed_Results/Fig4C_data.mat",'h1_vals','h2_vals','h3_vals','h4_vals','overall_mean1','overall_mean2','x_bins')
close all
%% FIGURE 4D: P(extinction) vs. relative abundance
clear all

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];
all_abuns = zeros(10000*500,1);
all_extinct = zeros(10000*500,1);
all_abuns2 = zeros(10000*500,1);
all_extinct2 = zeros(10000*500,1);
tot_all = 1;

R0=40;
fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(9)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
if size(deltas,1) ~= 10000
    size(deltas)
end
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    extinctions = d < -1e-3*std(h(1-h>0))/R0;
    extinctions(parent_inds(k)) = [];

    d = deltas{k,1};
    survivors = d > -1e-3*std(h(1-h>0))/R0;
    alphas = organisms{k,1}.enzymesInSpecies(survivors,:);
    fitnesses = organisms{k,1}.budget(survivors);
    alphas = alphas./sum(alphas,2);
    inv_harvest = 1./h;
    n_lsq = lsqnonneg(alphas', inv_harvest');
    n_lsq = n_lsq.*exp(-fitnesses);
    n_lsq = n_lsq/sum(n_lsq);

    extinctions = extinctions(1:end-1);
    n_lsq(parent_inds(k)) = [];

    all_abuns(tot_all:length(n_lsq)+tot_all-1) = n_lsq;
    all_extinct(tot_all:length(n_lsq)+tot_all-1) = extinctions;
    tot_all = tot_all + length(n_lsq);
end

edges = exp(-14:0.05:-2.5);
ext_probs = 0*edges(1:end-1);
ext_probs_err = 0*edges(1:end-1);
for i = 1:(length(edges)-1)
    slice = all_extinct(all_abuns > edges(i) & all_abuns < edges(i+1));
    ext_probs(i) = sum(slice)/length(slice);
    ext_probs_err(i) = sqrt(sum(slice))/length(slice);
end

save('Processed_Results/Fig4D_data.mat','all_abuns','all_extinct','edges','ext_probs','ext_probs_err')

%% FIGURE 4E: Fold change in use of resource targeted by mutation in extinct species
clear all

R0_vals = 10:10:100;
all_resource = -ones(10000*500,length(R0_vals),2);
all_extinct = -ones(10000*500,length(R0_vals),2);

for i = 1:length(R0_vals)
    tot=1;
    R0 = R0_vals(i);
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_" + string(R0_vals(i)) + "_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    clear organisms
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        extinctions = d < -1e-3*std(h(1-h>0))/R0;
        extinctions = extinctions(1:end-1);
        extinctions(parent_inds(k)) = [];
        enz_ind = find(organisms{k,2}.enzymesInSpecies(parent_inds(k),:)  ~= organisms{k,2}.enzymesInSpecies(end,:));
        if isempty(enz_ind)
            continue
        end

        enz_partial = organisms{k,2}.enzymesInSpecies(1:end-1,enz_ind);
        enz_partial(parent_inds(k)) = [];
        all_resource(tot:length(extinctions)+tot-1,i,1) = enz_partial;
        all_extinct(tot:length(extinctions)+tot-1,i,1) = extinctions;
        tot = tot + length(extinctions);

        if size(all_extinct, 1) > 10000*500
            size(all_extinct, 1)
        end
    end
end
for i = 1:length(R0_vals)
    tot=1;
    R0 = R0_vals(i);
    fname = "Raw_Results/First_Step_Sims/knockin/sstep_astar_0.8_phi_0.1_R0_" + string(R0_vals(i)) + "_R_200_dx_0_sigmaR2_0_muttype_0.mat";
    load(fname)
    if size(deltas,1) ~= 10000
        size(deltas)
    end
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        extinctions = d < -1e-3*std(h(1-h>0))/R0;
        extinctions = extinctions(1:end-1);
        extinctions(parent_inds(k)) = [];
        enz_ind = find(organisms{k,2}.enzymesInSpecies(parent_inds(k),:)  ~= organisms{k,2}.enzymesInSpecies(end,:));
        if isempty(enz_ind)
            continue
        end

        enz_partial = organisms{k,2}.enzymesInSpecies(1:end-1,enz_ind);
        enz_partial(parent_inds(k)) = [];
        all_resource(tot:length(extinctions)+tot-1,i,2) = enz_partial;
        all_extinct(tot:length(extinctions)+tot-1,i,2) = extinctions;
        tot = tot + length(extinctions);

        if size(all_extinct, 1) > 10000*500
            size(all_extinct, 1)
        end
    end
end

totals = zeros(4,10,2);
for i = 1:10
    for j = 1:2
        totals(1,i,j) = sum(all_extinct(:,i,j) == 0 & all_resource(:,i,j) == 0);
        totals(2,i,j) = sum(all_extinct(:,i,j) == 0 & all_resource(:,i,j) == 1);
        totals(3,i,j) = sum(all_extinct(:,i,j) == 1 & all_resource(:,i,j) == 0);
        totals(4,i,j) = sum(all_extinct(:,i,j) == 1 & all_resource(:,i,j) == 1);
    end
end

save('Processed_Results/Fig4E_data.mat','totals','R0_vals')

%% FIGURE 5A, right: Abundances at which different events occur
clear all

speciation_abun = zeros(1);
mutation_abun = zeros(1);
extinction_abun = zeros(1);
indices = [1 1 1];



for load_ind = 0:8

load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_" + string(load_ind) + ".mat")
species0 = 1:organisms_init.P;
last_ind = find(species_info(:,1),1);
speciation_counts = zeros(9000,1);
mutation_counts = zeros(9000,1);
extinction_counts = zeros(9000,1);


for i = 1:(length(all_abundances)-1)
    extinct_strains = setdiff(species_present{i},species_present{i+1});

    species1 = species_present{i};
    species2 = species_present{i+1};

    abuns1 = all_abundances{i};
    abuns2 = all_abundances{i+1};
    mut_ind = find(species1 == species_info(species2(end),2));
    for j = 1:length(species1)
        ind = find(species0 == species1(j));
        if j == length(species1)
            ind = length(abuns1);
        end
        if ~isempty(find(extinct_strains == species1(j),1)) && j ~= mut_ind
            %plot(i,abuns1(ind),'r.','MarkerSize',10)
            %hold on
            extinction_abun(indices(3)) = abuns1(ind);
            indices(3) = indices(3)+1;

            index = find_init_species(species1(mut_ind), species_info(:,2),last_ind-1);
            extinction_counts(index) = extinction_counts(index)+1;
        elseif j ~= mut_ind
            %plot([i i+1],[abuns1(ind) abuns2(j)],'Color',[0.8 0.8 0.8 0.5],'LineWidth',0.1)
            %hold on
        end
    end
    ind = find(species0 == species1(mut_ind));
    if mut_ind == length(species1)
        ind = length(abuns1);
    end
    if ~isempty(find(extinct_strains == species1(mut_ind),1))
        %plot([i i+1],[abuns1(ind) abuns2(length(species1))],'Color',[0.8 0.8 0.8 0.5],'LineWidth',0.1)
        %hold on
        %plot(i,abuns1(ind),'b.','MarkerSize',10)
        %hold on
        mutation_abun(indices(2)) = abuns1(ind);
        indices(2) = indices(2)+1;

        index = find_init_species(species1(mut_ind), species_info(:,2),last_ind-1);
        mutation_counts(index) = mutation_counts(index)+1;
    else
        %plot([i i+1],[abuns1(ind) abuns2(mut_ind)],'k','LineWidth',0.1)
        %hold on
        %plot([i i+1],[abuns1(ind) abuns2(length(species1))],'k','LineWidth',0.1)
        %plot(i,abuns1(ind),'g.','MarkerSize',10)
        speciation_abun(indices(1)) = abuns1(ind);
        indices(1) = indices(1)+1;

        index = find_init_species(species1(mut_ind), species_info(:,2),last_ind-1);
        speciation_counts(index) = speciation_counts(index)+1;
    end
    species0 = species1;
end
end

edges = -7:0.1:-1;

figure()
h1 = histogram(mutation_abun,'BinEdges',10.^(edges),'Normalization','probability');
figure()
h2 = histogram(extinction_abun,'BinEdges',10.^(edges),'Normalization','probability');
figure()
h3 = histogram(speciation_abun,'BinEdges',10.^(edges),'Normalization','probability');

edges = (edges(1:end-1) + edges(2:end))/2;
edges = 10.^edges;

save('Processed_Results/Fig5A_right_data.mat','h1','h2','h3','edges')
close all

%% FIGURE 5B, number of strains over time
clear all

N_species = zeros(9,1000);
N_strains = zeros(9,1000);
min_len = Inf;
max_len = 0;

for i = 0:8
    load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_" + string(i) + ".mat")
    min_len = min(min_len, length(all_abundances));
    max_len = max(max_len, length(all_abundances));
    for j = 1:length(all_abundances)
        N_strains(i+1,j) = length(species_present{j});
        last_ind = find(species_info(:,1),1);

        species_curr = species_present{j};
        original_parents = 0*species_curr;

        for k = 1:length(original_parents)
            original_parents(k) = find_init_species(species_curr(k), species_info(:,2),last_ind-1);
        end

        N_species(i+1,j) = length(unique(original_parents));
    end
end

for i = 1:9
    plot(N_strains(i,1:find(~N_strains(i,:),1)-1) - N_species(i,1:find(~N_species(i,:),1)-1),'Color',[0.5 0.5 1],'LineWidth',0.75)
    hold on
    p2=plot(mean(N_strains(:,1:min_len)) - mean(N_species(:,1:min_len)),'Color',[0 0 0.75],'LineWidth',2);
    plot(N_strains(i,1:find(~N_strains(i,:),1)-1),'Color',[1 0.5 0.5],'LineWidth',0.75)
    p1=plot(mean(N_strains(:,1:min_len)),'Color',[0.75 0 0],'LineWidth',2);
end

xlabel('No. mutations','FontSize',25)
ylabel('No. ecotypes','FontSize',25)
legend([p1 p2],{'Total','Closely related'})
set(gca,'FontSize',20)
xlim([0 max_len])

save('Processed_Results/Fig5B_data.mat','N_strains','N_species','min_len','max_len')

%% FIGURE 5C, surival time of coexisting relatives
clear all

coex_flag = -ones(9,1000);
coex_depths = -ones(9,1000);
h_stds = -ones(9,1000);
h_half_stds = -ones(9,1000);

running_coex = zeros(1,1000);
running_tot = zeros(1,1000);

for k = 0:8
    species0 = 1:4500;
    load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_" + string(k) + ".mat")
    for i = 1:(length(all_abundances)-1)
        extinct_strains = setdiff(species_present{i},species_present{i+1});
    
        species1 = species_present{i};
        species2 = species_present{i+1};
        abuns1 = all_abundances{i};
        abuns2 = all_abundances{i+1};
        mut_ind = find(species1 == species_info(species2(end),2));

        ind = find(species0 == species1(mut_ind));
        if mut_ind == length(species1)
            ind = length(abuns1);
        end
        if isempty(find(extinct_strains == species1(mut_ind),1))
            coex_flag(k+1,i)=1;
            coex_depths(k+1,i)=min(tree_depth(species1(mut_ind), i, species_present, species_info), tree_depth(species2(end), i, species_present, species_info));
        else
            coex_flag(k+1,i)=0;
        end
    
        species0 = species1;

        h = all_harvests{i};
        h_stds(k+1,i) = std(h);
        h_half_stds(k+1,i) = std(h(h<1));
    end

    for i = 1:length(par_abundances)
        running_coex(i) = running_coex(i) + mut_fitnesses(i) < (1 - 0.9)*par_abundances(i)*0.9*100/(20*20*0.9);
        running_tot(i) = running_tot(i) + 1;
    end
end

last_ind = zeros(10,1);
for i = 1:9
    last_ind(i) = find(coex_flag(i,:)==-1,1);
end

survival_prob = zeros(max(max(coex_depths)),1);
for i = 1:length(survival_prob)
    num = sum(sum(coex_depths >= i));
    den = 0;
    for j = 1:9
        den = den + sum(coex_flag(j,1:(last_ind(j)-i)) == 1);
    end
    survival_prob(i) = num/den;
end

save('Processed_Results/Fig5C_data.mat','survival_prob')

%% FIGURE 5D: DFE before and after evolution
clear all

R0 = 20;
R = 100;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.9;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot*sqrt(1 - alpha.*I_func(lmbda));

m = 0;
q = (R/R0)*(V_tot^2 - epsilon^2)/(1-R0/R);

x_vals = linspace(-0.01,0.01,1000);
S_inv_vals1 = normpdf(x_vals, m/R, sqrt(q/R)/R0);

fitnesses_ko1 = zeros(R^2*100,1);
ind_ko1 = 1;
fitnesses_ko2 = zeros(R^2*100,1);
ind_ko2 = 1;

for i = 0:8
    load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_" + string(i) + ".mat")

    h = all_harvests{1};
    h = length(h)*h/sum(h);
    enzymes = organisms_init.enzymesInSpecies;
    survivors = all_deltas{1} > -1e-3*std(h)/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    delta_alphas_ki = enzymes.*(1./(sum(enzymes,2)+1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    costs = organisms_init.budget;
    costs = costs(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(costs);
    overall_shift = -mean(log(overall_shift));
    costs = costs + overall_shift;

    for j = 1:size(enzymes,1)
        inv_fitness = sum(delta_alphas_ko(j,:).*h);
        inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
        fitnesses_ko1(ind_ko1:ind_ko1+length(inv_fitness)-1) = inv_fitness*exp(costs(j));
        ind_ko1 = ind_ko1+length(inv_fitness);
    end

    h = all_harvests{length(all_harvests)};
    h = length(h)*h/sum(h);
    enzymes = organisms_final.enzymesInSpecies;

    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    delta_alphas_ki = enzymes.*(1./(sum(enzymes,2)+1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    costs = organisms_final.budget;

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(costs);
    overall_shift = -mean(log(overall_shift));
    costs = costs + overall_shift;

    for j = 1:size(enzymes,1)
        inv_fitness = sum(delta_alphas_ko(j,:).*h);
        inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
        fitnesses_ko2(ind_ko2:ind_ko2+length(inv_fitness)-1) = inv_fitness*exp(costs(j));
        ind_ko2 = ind_ko2+length(inv_fitness);
    end
end

edges = (-0.5:0.002:0.5)/R0;
bins_ko_1 = zeros(length(edges)-1,1);
bins_ko_2 = zeros(length(edges)-1,1);

figure(5)
h1 = histogram(fitnesses_ko1(1:ind_ko1-1),'Normalization','pdf','BinEdges',edges);

figure(7)
h3 = histogram(fitnesses_ko2(1:ind_ko2-1),'Normalization','pdf','BinEdges',edges);

save('Processed_Results/Fig5D_data.mat','edges','h1','h3','x_vals','R0','S_inv_vals1')
close all

%% FIGURE 5E: P(coex) vs. evolutionary time
clear all

coex_flag = -ones(9,1000);
coex_depths = -ones(9,1000);
h_stds = -ones(9,1000);
h_half_stds = -ones(9,1000);

for k = 0:8
    species0 = 1:4500;
    load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_" + string(k) + ".mat")
    for i = 1:(length(all_abundances)-1)
        extinct_strains = setdiff(species_present{i},species_present{i+1});
    
        species1 = species_present{i};
        species2 = species_present{i+1};
        abuns1 = all_abundances{i};
        abuns2 = all_abundances{i+1};
        mut_ind = find(species1 == species_info(species2(end),2));

        ind = find(species0 == species1(mut_ind));
        if mut_ind == length(species1)
            ind = length(abuns1);
        end
        if isempty(find(extinct_strains == species1(mut_ind),1))
            coex_flag(k+1,i)=1;
            coex_depths(k+1,i)=min(tree_depth(species1(mut_ind), i, species_present, species_info), tree_depth(species2(end), i, species_present, species_info));
        else
            coex_flag(k+1,i)=0;
        end
    
        species0 = species1;

        h = all_harvests{i};
        h_stds(k+1,i) = std(h);
        h_half_stds(k+1,i) = std(h(h<1));
    end
end

coex = 0;
tot = 0;
R0=20;
fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:10000
    h = harvests{k,1};
    d = deltas{k,2};
    if d(end) < -1e-3*std(h(1-h>0))/R0
        continue
    end
    coex = coex + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
    tot = tot + 1;
end

cprob1 = coex/tot;

save('Processed_Results/Fig5E_data.mat','coex_flag','cprob1')

last_ind = zeros(10,1);
for i = 1:9
    last_ind(i) = find(coex_flag(i,:)==-1,1);
end

survival_prob = zeros(max(max(coex_depths)),1);
for i = 1:length(survival_prob)
    num = sum(sum(coex_depths >= i));
    den = 0;
    for j = 1:9
        den = den + sum(coex_flag(j,1:(last_ind(j)-i)) == 1);
    end
    survival_prob(i) = num/den;
end

function L = lmbda_func(lmbda, alpha, epsilon, R0, R)
    V_tot = (1-R0/R)*(2*(1 - alpha.*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
    L = V_tot.^2.*(1 - alpha*I_func(lmbda)) - epsilon.^2;
end

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end

function init_species = find_init_species(ind,parentage, original_N)
    if ind <= original_N
        init_species = ind;
    else
        init_species = find_init_species(parentage(ind),parentage,original_N);
    end
end

function depth = tree_depth(ind, gen, species_present, species_info)
    depth = 0;
    for i = (gen+1):length(species_present)
        if ~isempty(find(species_present{i}==ind,1))
            depth = depth + 1;
        else
            break
        end
    end

    daughter_species = find(species_info(:,1) > gen & species_info(:,2) == ind);
    for i = 1:length(daughter_species)
        depth = max(depth, species_info(daughter_species(i),1)-gen+tree_depth(daughter_species(i),species_info(daughter_species(i),1),species_present,species_info));
    end
end
