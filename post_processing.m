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
astar_actual = zeros(19,2,2);
astar_tots = zeros(19,2,2);
all_inds = ones(19,2,2);
astar_vals = 0.05:0.05:0.95;
phi_vals = [0.1 0.5];

R0 = 10;
for i = 1:19
    for i2 = 1:2
        fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_" + string(astar_vals(i)) + "_phi_" + string(phi_vals(i2)) + "_R0_10_R_200_sigmaR2_0.mat";
        load(fname)
        for i3 = 1:length(harvests)
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

            astar_actual(i,1,i2) = astar_actual(i,1,i2) + sum(survivors)/200;
            astar_tots(i,1,i2) = astar_tots(i,1,i2) + 1;

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
        for i3 = 1:length(harvests)
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

            astar_actual(i,2,i2) = astar_actual(i,2,i2) + sum(survivors)/200;
            astar_tots(i,2,i2) = astar_tots(i,2,i2) + 1;

            all_stds(i,2,i2,all_inds(i,2,i2)) = std(all_fitnesses(1:ind-1));
            all_inds(i,2,i2) = all_inds(i,2,i2)+1;
        end
    end
end

astar_actual = astar_actual./astar_tots;

sigma_invs = mean(all_stds,4);

save('Processed_Results/Fig2D_data','astar_actual','astar_vals','alpha_star','sigma_invs','sigma_inv1','sigma_inv2','sigma_inv3','sigma_inv4')

%% FIGURE 2E: Monoculture-community fitness correlation vs. niche saturation
clear all

astarvals = 0.1:0.1:0.9;
sigmaR2_vals = [0.007358 0.01213 0.02 0.03297 0.05437];
R0 = 40;
R = 200;
fitnesses_mono = zeros(9,5,R^2*1000);
fitnesses_comm = zeros(9,5,R^2*1000);
astar_actual = zeros(9,5);
astar_tots = zeros(9,5);
inds = ones(9,5);


for i = 1:9
    for i2 = 1:5
        fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_sigmaR2_" + string(sigmaR2_vals(i2)) + ".mat";
        load(fname)
        for i3 = 1:length(harvests)
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

            astar_actual(i,i2) = astar_actual(i,i2) + sum(survivors)/200;
            astar_tots(i,i2) = astar_tots(i,i2) + 1;
        end
    end
end
astar_actual = astar_actual./astar_tots;

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

save('Processed_Results/Fig2E_data','astarvals','astar_actual','alpha_star','all_corrs','corr_theory')

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

S_inv_all = zeros(100,2);
coex_all = false(100,2);

fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.3_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:size(organisms,1)
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
for k = 1:size(organisms,1)
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
for k = 1:size(organisms,1)
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
typical_shared_resources = typical_shared_resources(1:tot2-1);
num_shared_resources = num_shared_resources(1:tot_all);


tot_all = 0;
tot2=1;
R0 = 16;
fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(9)) + "_phi_0.1_R0_16_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:size(organisms,1)
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
        shared_resources = sum(organisms{k,2}.enzymesInSpecies(end,:).*other_enz,2);            
        typical_shared_resources2(tot2:tot2+length(shared_resources)-1) = shared_resources;
        tot2 = tot2 + length(shared_resources);
    end
end
typical_shared_resources2 = typical_shared_resources2(1:tot2-1);
num_shared_resources2 = num_shared_resources2(1:tot_all);


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
for k = 1:size(organisms,1)
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
    for k = 1:size(organisms,1)
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
    for k = 1:size(organisms,1)
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

%% FIGURE S1: Explicit community dynamics
clear all
addpath('Simulations')

seed = round(rand*100000000);

epsilon = 0.5;
alpha = 10;
sigma_R_sq = 0.1;
R = 100;
R0 = 20;

avg_supply = 5;

h_max = 4;
c_half = 1;

p.N = R;
p.R2 = sigma_R_sq*R;
p.sparsity = R0/R;
rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
while sum(p.capacity < 0) > 0
    fprintf('Negative resource supply generated, retrying.\n');
    p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
end
p.alpha = alpha;
p.seed = round(seed);
p.epsilon = epsilon/R0;

prop = generateOrganisms(p);
while min(sum(prop.enzymesInSpecies)) == 0
    fprintf('Unused resource, retrying.\n');
    p.seed = round(rand*100000000);
    prop = generateOrganisms(p);
end

prop.budget = prop.budget - max(prop.budget);
prop.capacity = p.capacity;


r = ( prop.enzymesInSpecies./sum(prop.enzymesInSpecies,2) ) .* exp(prop.budget);
[t,n,c] = log_community_dynamics(alpha*R,r,prop.capacity/mean(prop.capacity),avg_supply,h_max,c_half);

n(n<1e-10) = 0;

[h_eq, d_eq] = findEquilibrium(prop);

threshold = 1e-3*std(1-h_eq(1-h_eq > 0))/R0;
survivors = d_eq > -threshold;

alphas = prop.enzymesInSpecies(survivors,:);
fitnesses = prop.budget(survivors);
alphas = alphas./sum(alphas,2);
inv_harvest = prop.capacity'./(h_eq*mean(prop.capacity));
n_lsq = lsqnonneg(alphas', inv_harvest');
n_lsq = n_lsq.*exp(-fitnesses);
n_lsq = n_lsq/sum(n_lsq);

abundances = zeros(length(d_eq),1);
abundances(survivors) = n_lsq;

h1 = harvest(c(end,:),h_max,c_half);
h1 = h1/sum(h1);
h2 = h_eq/sum(h_eq);

[h_dynamics, inds] = sort(h1/mean(h1));
h2 = h2(inds);
[n_dynamics, inds] = sort(n(end,:)/sum(n(end,:)));
live_flag = n_dynamics > 0 | abundances(inds)' > 0;

abundances = abundances(inds(live_flag));
n_dynamics = n_dynamics(live_flag);

h2_1 = h2;
h_dynamics_1 = h_dynamics;
abundances_1 = abundances;
n_dynamics_1 = n_dynamics;
t1 = t;
n1 = n;
c1 = c;

seed = round(rand*100000000);

epsilon = 0.5;
alpha = 10;
sigma_R_sq = 0.1;
R = 100;
R0 = 20;

avg_supply = 50;

h_max = 4;
c_half = 1;

p.N = R;
p.R2 = sigma_R_sq*R;
p.sparsity = R0/R;
rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
while sum(p.capacity < 0) > 0
    fprintf('Negative resource supply generated, retrying.\n');
    p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N).*normrnd(0,1,p.N,1));
end
p.alpha = alpha;
p.seed = round(seed);
p.epsilon = epsilon/R0;

prop = generateOrganisms(p);
while min(sum(prop.enzymesInSpecies)) == 0
    fprintf('Unused resource, retrying.\n');
    p.seed = round(rand*100000000);
    prop = generateOrganisms(p);
end

prop.budget = prop.budget - max(prop.budget);
prop.capacity = p.capacity;


r = ( prop.enzymesInSpecies./sum(prop.enzymesInSpecies,2) ) .* exp(prop.budget);
[t,n,c] = log_community_dynamics(alpha*R,r,prop.capacity/mean(prop.capacity),avg_supply,h_max,c_half);

n(n<1e-10) = 0;

[h_eq, d_eq] = findEquilibrium(prop);

threshold = 1e-3*std(1-h_eq(1-h_eq > 0))/R0;
survivors = d_eq > -threshold;

alphas = prop.enzymesInSpecies(survivors,:);
fitnesses = prop.budget(survivors);
alphas = alphas./sum(alphas,2);
inv_harvest = prop.capacity'./(h_eq*mean(prop.capacity));
n_lsq = lsqnonneg(alphas', inv_harvest');
n_lsq = n_lsq.*exp(-fitnesses);
n_lsq = n_lsq/sum(n_lsq);

abundances = zeros(length(d_eq),1);
abundances(survivors) = n_lsq;

h1 = harvest(c(end,:),h_max,c_half);
h1 = h1/sum(h1);
h2 = h_eq/sum(h_eq);

[h_dynamics, inds] = sort(h1/mean(h1));
h2 = h2(inds);

[n_dynamics, inds] = sort(n(end,:)/sum(n(end,:)));
live_flag = n_dynamics > 0 | abundances(inds)' > 0;

abundances = abundances(inds(live_flag));
n_dynamics = n_dynamics(live_flag);


h2_2 = h2;
h_dynamics_2 = h_dynamics;
abundances_2 = abundances;
n_dynamics_2 = n_dynamics;

save('Processed_Results/FigS1_data.mat','t1','n1','c1','abundances_1','n_dynamics_1','c_half','h_max','h2_1','h_dynamics_1','abundances_2','n_dynamics_2','h2_2','h_dynamics_2')


%% FIGURE S2: DFE for single organism
clear all
R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = 0.8;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;

x_vals = linspace(-0.035,0.035,1000);
S_inv_vals3 = normpdf(x_vals, 0, sigma_inv);

fname = "Raw_Results/Assembly_Sims/astar_phi_params/comm_astar_0.8_phi_0.1_R0_40_R_200_sigmaR2_0.mat";
load(fname)
for i = 1:1
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
    for j = 1:1
        inv_fitness = sum(delta_alphas_ko(j,:).*h);
        inv_fitness = inv_fitness - (delta_alphas_ko(j,find(enzymes(j,:),1)) + 1/(sum(enzymes(j,:))))*h(enzymes(j,:));
        fitnesses_ko3 = inv_fitness*exp(budgets(j));

        inv_fitness = sum(delta_alphas_ki(j,:).*h);
        inv_fitness = inv_fitness + (1/(1+sum(enzymes(j,:))))*h(~enzymes(j,:));
        fitnesses_ki3 = inv_fitness*exp(budgets(j));
    end
end

edges2 = (-0.0026:0.0004:0.0026);

figure()
h5 = histogram(fitnesses_ko3,'Normalization','pdf','BinEdges',edges2);

figure()
h6 = histogram(fitnesses_ki3,'Normalization','pdf','BinEdges',edges2);
save('Processed_Results/FigS2_data.mat','edges2','h5','h6','x_vals','S_inv_vals3')
close all

%% Figure S3A and S3B: Resource availabilities and P(coex) for Dirichlet uptake rates
clear all

coex = zeros(14,1);
tot = zeros(14,1);
std_h = zeros(14,1);
harvests1 = zeros(200*10000,1);
harvests2 = zeros(200*10000,1);
h_ind1 = 1;
h_ind2 = 1;


astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];

astar_actual = zeros(14,1);

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_alttype_2_gamma_1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual(i) = astar_actual(i) + (length(d)-1)/200;
        coex(i) = coex(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        std_h(i) = std_h(i) + std(h);
        tot(i) = tot(i) + 1;

        if i == 5
            harvests1(h_ind1:h_ind1+199) = h;
            h_ind1 = h_ind1+200;
        elseif i == 8
            harvests2(h_ind2:h_ind2+199) = h;
            h_ind2 = h_ind2+200;
        end
    end
end

R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.8;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot*sqrt(1 - alpha.*I_func(lmbda));

q = (R/R0)*(V_tot^2 - epsilon^2)/(1-R0/R);

x_vals = linspace(0.5,1.5,1000);
S_inv_vals2 = normpdf(x_vals, 1, sqrt(q/R));


alpha_star = 0.5;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot*sqrt(1 - alpha.*I_func(lmbda));

q = (R/R0)*(V_tot^2 - epsilon^2)/(1-R0/R);
S_inv_vals1 = normpdf(x_vals, 1, sqrt(q/R));

edges = 0.5:0.001:1.5;

figure()
h1=histogram(harvests1,'BinEdges',edges,'Normalization','pdf');
figure()
h2=histogram(harvests2,'BinEdges',edges,'Normalization','pdf');

h1Vals = h1.Values;
h2Vals = h2.Values;
save('Processed_Results/FigS3A_data','edges','h1Vals','h2Vals','x_vals','S_inv_vals1','S_inv_vals2')
close all

R0 = 40;
R = 200;

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.001:0.001:0.999;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

save('Processed_Results/FigS3B_data','alpha_star_vals','astar_actual', 'cprob', 'astarvals','coex','tot')

%% FIGURE S3C: P(coex) for Gaussian resource uptakes
clear all

coex1 = zeros(8,1);
tot1 = zeros(8,1);
astar_actual1 = zeros(8,1);

coex2 = zeros(8,1);
tot2 = zeros(8,1);
astar_actual2 = zeros(8,1);

coex3 = zeros(8,1);
tot3 = zeros(8,1);
astar_actual3 = zeros(8,1);

coex4 = zeros(8,1);
tot4 = zeros(8,1);
astar_actual4 = zeros(8,1);

var_alphas = zeros(8,4);
tot_var = zeros(8,4);

all_harvests = zeros(200*1000,8,4);
h_inds = ones(8,4);

alphavals = [0.3 0.6 1.2 2.4 4.8 9.6 19.2 38.4];

R0=10*200;
R = 200;
for i = 1:8
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_gauss_uptakes_alpha_" + string(alphavals(i)) + "_eps_1_scale_10_R_200_mutsize_0.1.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/R;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot1(i) = tot1(i) + 1;
    end
    all_harvests(h_inds(i,1):h_inds(i,1)+length(h)-1,i,1) = h;
    h_inds(i,1) = h_inds(i,1) + length(h);
end

R0=100*200;
R = 200;
for i = 1:8
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_gauss_uptakes_alpha_" + string(alphavals(i)) + "_eps_1_scale_100_R_200_mutsize_0.1.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/R;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot2(i) = tot2(i) + 1;
    end
end

R0=10*200;
R = 200;
for i = 1:8
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_gauss_uptakes_alpha_" + string(alphavals(i)) + "_eps_1_scale_10_R_200_mutsize_0.05.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual3(i) = astar_actual3(i) + (length(d)-1)/R;
        coex3(i) = coex3(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot3(i) = tot3(i) + 1;
    end
end

R0=100*200;
R = 200;
for i = 1:8
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_gauss_uptakes_alpha_" + string(alphavals(i)) + "_eps_1_scale_100_R_200_mutsize_0.05.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual4(i) = astar_actual4(i) + (length(d)-1)/R;
        coex4(i) = coex4(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot4(i) = tot4(i) + 1;
    end
end

N_val = 200;
B = 100;
mut_size = 0.1;

R0 = N_val/(1+1/B);

phi = (astar_actual2./tot2)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.2:0.01:0.8;

lmbda = interp1(astar_actual2./tot2,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val/(R0*(1-R0/N_val));

sinv_factor = sqrt(B*q/mut_size)/R0^2;
Csq = (V_tot./(R0^2*(1-R0/N_val)*sinv_factor)).^2;



L = lmbda;
cprob = (-2.*(Csq + Csq.^2) + lmbda.*sqrt(2.*pi).*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)./2.).*erfc(lmbda./sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2.*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)./2.).*lmbda.*sqrt(2.*pi).*erfc(lmbda./sqrt(2))));


mut_size = 0.05;

R0 = N_val/(1+1/B);

phi = (astar_actual4./tot4)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.2:0.01:0.8;

lmbda = interp1(astar_actual4./tot4,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val/(R0*(1-R0/N_val));

sinv_factor = sqrt(B*q/mut_size)/R0^2;
Csq = (V_tot./(R0^2*(1-R0/N_val)*sinv_factor)).^2;



L = lmbda;
cprob2 = (-2.*(Csq + Csq.^2) + lmbda.*sqrt(2.*pi).*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)./2.).*erfc(lmbda./sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2.*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)./2.).*lmbda.*sqrt(2.*pi).*erfc(lmbda./sqrt(2))));

save('Processed_Results/FigS3C_data','alpha_star_vals','cprob','cprob2','astar_actual1','astar_actual2','astar_actual3','astar_actual4','coex1','coex2','coex3','coex4','tot1','tot2','tot3','tot4')

%% FIGURE S3D: Specialist community assembly
clear all
addpath('Simulations')
seed = round(rand*10000000);
R = 50;
epsilon = 5;
gamma = 0.1;

R2_val = 1*R;

p.N = R;
p.R2 = R2_val;
p.sparsity = 1/R;
p.capacity = 100*exprnd(1,p.N,1);
p.alpha = 1;
p.seed = seed;
p.epsilon = epsilon;
prop = generateOrganismsSpecialists(p);
prop.capacity = p.capacity/mean(p.capacity);
true_deltas = {};

[h, delta] = findEquilibrium(prop);
budgets = prop.budget;

invasion_flag = zeros(R);
coex_flag = zeros(R);
all_deltas = zeros(R*R,1);
ind = 1;

for i = 1:R
    for j = 1:R
        if budgets(i) > budgets(j) && i ~= j
            invasion_flag(i,j) = 1;
            prop2 = prop;
            prop2.enzymesInSpecies = [eye(R); zeros(1,R)];
            prop2.enzymesInSpecies(end,i) = 1 - gamma;
            prop2.enzymesInSpecies(end,j) = gamma;
            prop2.budget = [budgets; budgets(i)];

            [h2, delta2] = findEquilibrium(prop2);
            if delta2(i) > -1e-12
                coex_flag(i,j) = 1;
            end
            if delta2(end) < -1e-10
                delta2(end)
            end
            all_deltas(ind) = delta2(i);
            ind = ind+1;
        end
    end
end
all_deltas = all_deltas(1:ind-1);
histogram(log10(-all_deltas))

invasion_flag = invasion_flag(:);
coex_flag = coex_flag(:);

Ksource_vals = repmat(prop.capacity,1,R);
Krecip_vals = repmat(prop.capacity',R,1);
Ksource_vals = Ksource_vals(:);
Krecip_vals = Krecip_vals(:);

save('Processed_Results/FigS3D_data','Ksource_vals','Krecip_vals','coex_flag','invasion_flag','gamma','prop','R')

%% FIGURE S4: "Simultaneous assembly" simulations without extinction
clear all

coex1 = zeros(14,1);
tot1 = zeros(14,1);
coex2 = zeros(14,1);
tot2 = zeros(14,1);

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];

astar_actual1 = zeros(14,1);
astar_actual2 = zeros(14,1);

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/200;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot1(i) = tot1(i) + 1;
    end
end

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/simul_assembly/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            d(end)
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (sum(deltas{k,1} > -1e-3*std(h(1-h>0))/R0))/200;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
R = 200;

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.001:0.001:0.999;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));


save('Processed_Results/FigS4_data','alpha_star_vals','astar_actual1','astar_actual2', 'cprob', 'astarvals','coex1','coex2','tot1','tot2')

%% FIGURE S5, left: P(coex) vs niche saturation for knockout and global-effect mutations
clear all

coex1 = zeros(14,1);
tot1 = zeros(14,1);
coex2 = zeros(14,1);
tot2 = zeros(14,1);

astarvals = 0.1:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];

astar_actual1 = zeros(14,1);
astar_actual2 = zeros(14,1);

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_alttype_2_gamma_1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/200;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot1(i) = tot1(i) + 1;
    end
end


R0=40;
 for i = 1:14
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_alttype_3_gamma_1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            d(end)
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/200;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
R = 200;

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.001:0.001:0.999;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));


save('Processed_Results/FigS5_left_data','alpha_star_vals','astar_actual1','astar_actual2', 'cprob', 'astarvals','coex1','coex2','tot1','tot2')

%% FIGURE S5, right: P(coex) vs gamma for knockout and global-effect mutations
clear all

coex1 = zeros(10,1);
tot1 = zeros(10,1);
coex2 = zeros(10,1);
tot2 = zeros(10,1);

gammavals = 0.1:0.1:1;

R0=40;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_astar_0.8_phi_0.1_R0_40_R_200_alttype_2_gamma_" + string(gammavals(i)) +".mat";
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

R0=40;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/dirichlet_uptakes/sstep_astar_0.8_phi_0.1_R0_40_R_200_alttype_3_gamma_" + string(gammavals(i)) +".mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            d(end)
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
alpha_star_vals = 0.8;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

save('Processed_Results/FigS5_right_data','gammavals', 'cprob','coex1','coex2','tot1','tot2')

%% FIGURE S6: Invasion fitness vs. parent rel. abun. scatterplot
clear all
rel_abuns = [];
s_invs = [];
coex_flags = [];
R0 = 40;

fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
load(fname)
for k = 1:10000
    h1 = harvests{k,1};
    threshold = 1e-3*std(1-h1(1-h1 > 0))/R0;
    survivors = deltas{k,1} > - threshold;
    inv_harvest = 1./h1;
    alphas = organisms{k,1}.enzymesInSpecies(survivors,:);
    fitnesses = organisms{k,1}.budget(survivors);
    alphas = alphas./sum(alphas,2);
    n_lsq = lsqnonneg(alphas',inv_harvest');
    n_lsq = n_lsq.*exp(-fitnesses);
    n_lsq = n_lsq/sum(n_lsq);

    rel_abuns(end+1) = n_lsq(parent_inds(k));

    enz_ind = find(organisms{k,2}.enzymesInSpecies(parent_inds(k),:) ~= organisms{k,2}.enzymesInSpecies(end,:));
    alphas_old = organisms{k,2}.enzymesInSpecies(parent_inds(k),:);
    alphas_new = organisms{k,2}.enzymesInSpecies(end,:);
    alphas_old = alphas_old/sum(alphas_old);
    alphas_new = alphas_new/sum(alphas_new);

    h = length(h1)*h1/sum(h1);
    enzymes = organisms{k,1}.enzymesInSpecies;
    survivors = deltas{k,1} > -1e-3*std(1-h1(1-h1 > 0))/R0;
    enzymes = enzymes(survivors,:);
    delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
    alphas = enzymes.*(1./sum(enzymes,2));
    budgets = organisms{k,1}.budget;
    budgets = budgets(survivors);

    overall_shift = sum(alphas.*h,2);
    overall_shift = overall_shift.*exp(budgets);
    overall_shift = -mean(log(overall_shift));
    budgets = budgets + overall_shift;

    s_invs(end+1) = sum(h1.*(alphas_new-alphas_old))*exp(budgets(parent_inds(k)));

    deltas2 = deltas{k,2};
    coex_flags(end+1) = deltas2(parent_inds(k)) > -threshold;
end

coex_flags = ~~coex_flags;
save('Processed_Results/FigS6_data','s_invs','rel_abuns','coex_flags');

%% FIGURE S7A P(Coex) for small deviations from uniform resource supply
clear all

coex1 = zeros(9,1);
tot1 = zeros(9,1);

sigmaR2_vals = 0.01:0.01:0.09;

astar_actual1 = zeros(9,1);

R0=40;
for i = 1:9
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_0.8_phi_0.1_R0_40_R_200_dx_0_sigmaR2_" + string(sigmaR2_vals(i)) + "_muttype_-1.mat";

    load(fname)

    for k = 1:length(harvests)
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/200;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0))/R0);
        tot1(i) = tot1(i) + 1;
    end
end

R0 = 40;
N_val = 200;

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.8;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
sigmaR2_vals = 0:0.001:0.095;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val/(R0*(1-R0/N_val)) + N_val*sigmaR2_vals.*(1-alpha_star_vals).^2;
sigma_inv = sqrt(q/N_val)/R0;
Csq = (V_tot./(R0^2*(1-R0/N_val)*sigma_inv)).^2;

cprob = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

save('Processed_Results/FigS7A_data','coex1','tot1','sigmaR2_vals','cprob');

%% FIGURE S7B, S7C P(coex) and DFE for exponential resource supply
clear all

coex1 = zeros(10,1);
tot1 = zeros(10,1);
astar_actual1 = zeros(10,1);

coex2 = zeros(10,1);
tot2 = zeros(10,1);
astar_actual2 = zeros(10,1);

coex3 = zeros(10,1);
tot3 = zeros(10,1);
astar_actual3 = zeros(10,1);

alphavals = [0.3 0.6 1.2 2.4 4.8 9.6 19.2 38.4 76.8 153.6];

fitnesses_ko = zeros(200^2*1000,10,5);
ind_ko = ones(10,5);

R0=30;
R = 150;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/exponential_supply/sstep_wide_resource_var_alpha_" + string(alphavals(i)) + "_eps_0.1_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/R;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot1(i) = tot1(i) + 1;

        if k <= length(organisms)
            enzymes = organisms{k,1}.enzymesInSpecies;
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
    
            for ind = 1:size(enzymes,1)
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,1):ind_ko(i,1)+length(inv_fitness)-1,i,1) = inv_fitness*exp(budgets(ind));
                ind_ko(i,1) = ind_ko(i,1)+length(inv_fitness);
            end
        end
    end
end

for i = 1:10
    fname = "Raw_Results/First_Step_Sims/exponential_supply/sstep_wide_resource_var_alpha_" + string(alphavals(i)) + "_eps_1_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/R;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot2(i) = tot2(i) + 1;

        if k <= length(organisms)
            enzymes = organisms{k,1}.enzymesInSpecies;
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
    
            for ind = 1:size(enzymes,1)
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,2):ind_ko(i,2)+length(inv_fitness)-1,i,2) = inv_fitness*exp(budgets(ind));
                ind_ko(i,2) = ind_ko(i,2)+length(inv_fitness);
            end
        end
    end
end

for i = 1:10
    fname = "Raw_Results/First_Step_Sims/exponential_supply/sstep_wide_resource_var_alpha_" + string(alphavals(i)) + "_eps_5_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual3(i) = astar_actual3(i) + (length(d)-1)/R;
        coex3(i) = coex3(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot3(i) = tot3(i) + 1;

        if k <= length(organisms)
            enzymes = organisms{k,1}.enzymesInSpecies;
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;
    
            for ind = 1:size(enzymes,1) 
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,3):ind_ko(i,3)+length(inv_fitness)-1,i,3) = inv_fitness*exp(budgets(ind));
                ind_ko(i,3) = ind_ko(i,3)+length(inv_fitness);
            end
        end
    end
end

fitnesses_ko = fitnesses_ko(1:max(max(ind_ko))-1,:,:);

R0 = 30;
N_val = 150;

phi = (astar_actual2./tot2)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.01:0.01:0.8;
lmbda = interp1(astar_actual2./tot2,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;
sigmaR2_vals = 1;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val/(R0*(1-R0/N_val)) + N_val*sigmaR2_vals.*(1-alpha_star_vals).^2;
sigma_inv = sqrt(q/N_val)/R0;
Csq = (V_tot./(R0^2*(1-R0/N_val)*sigma_inv)).^2;

cprob = (-2.*(Csq + Csq.^2) + lmbda.*sqrt(2.*pi).*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)./2.).*erfc(lmbda./sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2.*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)./2.).*lmbda.*sqrt(2.*pi).*erfc(lmbda./sqrt(2))));

save('Processed_Results/FigS7B_data.mat','astar_actual1','astar_actual2','astar_actual3','tot1','tot2','tot3','coex1','coex2','coex3','alpha_star_vals','cprob')

i1 = 4;
i2 = 10;

R0 = 30;
R = 150;
phi = astar_actual2(i1)./tot2(i1) / alphavals(i1);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual2(i1)./tot2(i1);
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
x_vals = linspace(-0.2,0.2,1000);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
norm_pdf = normpdf(x_vals,0,sigma_inv);
exp_pdf = (R0/(1-alpha_star))*exp(-(1+R0*x_vals/(1-alpha_star))).*(x_vals > -(1-alpha_star)/R0);

S_inv_vals1 = conv(norm_pdf,exp_pdf,'same');
S_inv_vals1 = flip((S_inv_vals1/sum(S_inv_vals1))/(x_vals(2)-x_vals(1)));

q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R)) + R*(1-alpha_star)^2;
S_inv_vals1B = normpdf(x_vals,0,sqrt(q/R)/R0);

R0 = 30;
R = 150;
phi = astar_actual2(i2)./tot2(i2) / alphavals(i2);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual2(i2)./tot2(i2);
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
norm_pdf = normpdf(x_vals,0,sigma_inv);
exp_pdf = (R0/(1-alpha_star))*exp(-(1+R0*x_vals/(1-alpha_star))).*(x_vals > -(1-alpha_star)/R0);

S_inv_vals2 = conv(norm_pdf,exp_pdf,'same');
S_inv_vals2 = flip((S_inv_vals2/sum(S_inv_vals2))/(x_vals(2)-x_vals(1)));

q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R)) + R*(1-alpha_star)^2;
S_inv_vals2B = normpdf(x_vals,0,sqrt(q/R)/R0);


figure()
h1 = histogram(fitnesses_ko(1:ind_ko(i1,2)-1,i1,2),'Normalization','pdf');
figure()
h2 = histogram(fitnesses_ko(1:ind_ko(i2,2)-1,i2,2),'Normalization','pdf');

%plot(x_vals,S_inv_vals1,'k--','LineWidth',1.25)
%plot(x_vals,S_inv_vals2,'r--','LineWidth',1.25)

save('Processed_Results/FigS7C_data.mat','x_vals','S_inv_vals1B','S_inv_vals2B','h1','h2')
close all

%% FIGURE S8A: DFE for Tikhonov & Monasson consumer resource model
clear all
R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.8;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
epsilon = V_tot*sqrt(1 - alpha.*I_func(lmbda));

m = (R/R0)*V_tot*lmbda;
q = (R/R0)*(V_tot^2 - epsilon^2)/(1-R0/R);

x_vals = linspace(-1,1,1000);
S_inv_vals1 = normpdf(x_vals, m/R, sqrt(q/R));
S_inv_vals2 = normpdf(x_vals, -m/R, sqrt(q/R));

fitnesses_ko = zeros(200*200*1000,1);
ind_ko = 1;

fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_astar_0.8_phi_0.1_R0_40_R_200_alttype_0.mat";
load(fname)

for i = 1:1000
    h = harvests{i};
    enzymes = organisms{i}.enzymesInSpecies;
    survivors = deltas{i} > -1e-3*std(h);
    enzymes = enzymes(survivors,:);
    for j = 1:size(enzymes,1)
        inv_fitness = 1 - h(~~enzymes(j,:));
        fitnesses_ko(ind_ko:ind_ko+length(inv_fitness)-1) = inv_fitness;
        ind_ko = ind_ko + length(inv_fitness);
    end
    fitnesses_ko((i-1)*200+1:i*200) = 1-harvests{i}';
end
fitnesses_ko = fitnesses_ko(1:ind_ko-1);

edges = (-0.2:0.002:0.2);
bins_ko = zeros(length(edges)-1,1);
bins_ki = zeros(length(edges)-1,1);

figure(5)
h1 = histogram(fitnesses_ko,'Normalization','pdf','BinEdges',edges);

figure(6)
h2 = histogram(-fitnesses_ko,'Normalization','pdf','BinEdges',edges);

save('Processed_Results/FigS8A_data.mat','edges','h1','h2','x_vals','R0','S_inv_vals1')
close all

%% FIGURE S8B: P(coex) for knockins and knockouts in Tikhonov model
clear all

coex1 = zeros(13,1);
tot1 = zeros(13,1);

coex2 = zeros(13,1);
tot2 = zeros(13,1);


astarvals = 0.2:0.1:0.9;
astarvals = [astarvals 0.92 0.94 0.96 0.98 0.99];

astar_actual1 = zeros(13,1);
astar_actual2 = zeros(13,1);

R0=40;
for i = 1:13
    fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_alttype_0.mat";
    load(fname)
    for k = 1:size(deltas,1)
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/200;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0)));
        tot1(i) = tot1(i) + 1;
    end
end

for i = 1:13
    fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_alttype_1.mat";
    load(fname)
    for k = 1:size(deltas,1)
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/200;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h(1-h>0)));
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.01:0.001:0.99;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/R)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
m = (R/R0)*V_tot*lmbda;
q = V_tot.^2.*alpha_vals.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
dx_vals = m/(R)/R0;

Csq = (V_tot./(R0^2*(1-R0/R)*sigma_inv)).^2;

cprob1 = (-2*(Csq + Csq.^2) + lmbda*sqrt(2*pi)*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)/2.).*erfc(lmbda/sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)/2.).*lmbda*sqrt(2*pi).*erfc(lmbda/sqrt(2))));

fun = @ (C,dx,sigma) (sigma.*(2.*(-power(exp(1),-power(C - dx,2)./(2..*power(sigma,2))) + power(exp(1),-power(dx,2)./(2.*power(sigma,2)))).*sigma + dx.*sqrt(2.*pi).*erf(dx./(sqrt(2).*sigma)) -... 
       dx.*sqrt(2.*pi).*erf((-C + dx)./(sqrt(2).*sigma))))./...
   (2.*(power(sigma,2)./power(exp(1),power(dx,2)./(2.*power(sigma,2))) + dx.*sqrt(pi./2).*sigma.*(1 + erf(dx./(sqrt(2).*sigma)))));

fun2 = @ (gma, R_0, R, Vtot, lma, n, dx, sigma) (-n - lma).*exp(-n.^2./2)./(exp(-lma.^2./2)-lma.*sqrt(pi./2).*erfc(lma./sqrt(2))) .* fun(gma.^2.*Vtot.*(-n-lma)./(R_0.^2.*(1-R_0./R)), dx, sigma);

pcoex_vals1 = 0*alpha_star_vals;
pcoex_vals2 = 0*alpha_star_vals;
pcoex_vals3 = 0*alpha_star_vals;
for i = 1:length(alpha_star_vals)
    pcoex_vals1(i) = integral(@(x) fun2(1,R0,R,V_tot(i),lmbda,x,dx_vals(i),sigma_inv(i)),-10,-lmbda);
    pcoex_vals2(i) = integral(@(x) fun2(1,R0,R,V_tot(i),lmbda,x,-dx_vals(i),sigma_inv(i)),-10,-lmbda);
    pcoex_vals3(i) = integral(@(x) fun2(1,R0,R,V_tot(i),lmbda,x,0*dx_vals(i),sigma_inv(i)),-10,-lmbda);
end

save('Processed_Results/FigS8B_data.mat','alpha_star_vals','astar_actual1','astar_actual2','pcoex_vals1','pcoex_vals2','pcoex_vals3','coex1','tot1','coex2','tot2')

%% FIGURE S8C, S8D: DFE and P(Coex) for communities without metabolic tradeoffs
clear all

coex1 = zeros(10,1);
tot1 = zeros(10,1);
astar_actual1 = zeros(10,1);

coex2 = zeros(10,1);
tot2 = zeros(10,1);
astar_actual2 = zeros(10,1);

coex3 = zeros(10,1);
tot3 = zeros(10,1);
astar_actual3 = zeros(10,1);

alphavals = [0.3 0.6 1.2 2.4 4.8 9.6 19.2 38.4 76.8 153.6];

fitnesses_ko = zeros(200^2*1000,10,3);

ind_ko = ones(10,3);

R0_actual = zeros(200*1000,10,3);
ind_R0 = ones(10,3);


R0 = 30;
R = 150;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_no_tradeoffs_alpha_" + string(alphavals(i)) + "_eps_0.1_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/R;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot1(i) = tot1(i) + 1;

        if k <= length(organisms)

        enzymes = logical(organisms{k,1}.enzymesInSpecies);
        survivors = deltas{k,1} > -1e-3*std(h)/R0;
        enzymes = enzymes(survivors,:);
        delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
        alphas = enzymes.*(1./sum(enzymes,2));
        budgets = organisms{k,1}.budget;
        budgets = budgets(survivors);
    
        overall_shift = sum(alphas.*h,2);
        overall_shift = overall_shift.*exp(budgets);
        overall_shift = -mean(log(overall_shift));
        budgets = budgets + overall_shift;

        for ind = 1:size(enzymes,1)
            inv_fitness = sum(delta_alphas_ko(ind,:).*h);
            inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
            fitnesses_ko(ind_ko(i,1):ind_ko(i,1)+length(inv_fitness)-1,i,1) = inv_fitness*exp(budgets(ind));
            ind_ko(i,1) = ind_ko(i,1)+length(inv_fitness);

            R0_actual(ind_R0(i,1),i,1) = sum(enzymes(ind,:));
            ind_R0(i,1) = ind_R0(i,1) + 1;
        end
        end
    end
end

for i = 1:10
    fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_no_tradeoffs_alpha_" + string(alphavals(i)) + "_eps_1_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/R;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot2(i) = tot2(i) + 1;

        if k <= length(organisms)

        enzymes = logical(organisms{k,1}.enzymesInSpecies);
        survivors = deltas{k,1} > -1e-3*std(h)/R0;
        enzymes = enzymes(survivors,:);
        delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
        alphas = enzymes.*(1./sum(enzymes,2));
        budgets = organisms{k,1}.budget;
        budgets = budgets(survivors);
    
        overall_shift = sum(alphas.*h,2);
        overall_shift = overall_shift.*exp(budgets);
        overall_shift = -mean(log(overall_shift));
        budgets = budgets + overall_shift;

        for ind = 1:size(enzymes,1)
            inv_fitness = sum(delta_alphas_ko(ind,:).*h);
            inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
            fitnesses_ko(ind_ko(i,2):ind_ko(i,2)+length(inv_fitness)-1,i,2) = inv_fitness*exp(budgets(ind));
            ind_ko(i,2) = ind_ko(i,2)+length(inv_fitness);

            R0_actual(ind_R0(i,2),i,2) = sum(enzymes(ind,:));
            ind_R0(i,2) = ind_R0(i,2) + 1;
        end
        end
    end
end

for i = 1:10
    fname = "Raw_Results/First_Step_Sims/metabolic_tradeoffs/sstep_no_tradeoffs_alpha_" + string(alphavals(i)) + "_eps_5_R0_30_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual3(i) = astar_actual3(i) + (length(d)-1)/R;
        coex3(i) = coex3(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot3(i) = tot3(i) + 1;

        if k <= length(organisms)

        enzymes = logical(organisms{k,1}.enzymesInSpecies);
        survivors = deltas{k,1} > -1e-3*std(h)/R0;
        enzymes = enzymes(survivors,:);
        delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
        alphas = enzymes.*(1./sum(enzymes,2));
        budgets = organisms{k,1}.budget;
        budgets = budgets(survivors);
    
        overall_shift = sum(alphas.*h,2);
        overall_shift = overall_shift.*exp(budgets);
        overall_shift = -mean(log(overall_shift));
        budgets = budgets + overall_shift;

        for ind = 1:size(enzymes,1)
            inv_fitness = sum(delta_alphas_ko(ind,:).*h);
            inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
            fitnesses_ko(ind_ko(i,3):ind_ko(i,3)+length(inv_fitness)-1,i,3) = inv_fitness*exp(budgets(ind));
            ind_ko(i,3) = ind_ko(i,3)+length(inv_fitness);

            R0_actual(ind_R0(i,3),i,3) = sum(enzymes(ind,:));
            ind_R0(i,3) = ind_R0(i,3) + 1;
        end
        end
    end
end

fitnesses_ko = fitnesses_ko(1:max(max(ind_ko))-1,:,:);
R0_actual = R0_actual(1:max(max(ind_R0))-1,:,:);

R0 = mean(R0_actual(1:ind_R0(:,3)-1,:,3));
N_val = 150;

phi = (astar_actual2./tot2)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.1:0.01:0.5;
lmbda = interp1(astar_actual1./tot1,lmbda,alpha_star_vals,'linear','extrap');
R0 = interp1(astar_actual1./tot1,R0,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0./N_val).*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
cost = 0;

q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val./(R0.*(1-R0./N_val));
sigma_inv = sqrt(q/N_val)./R0;
Csq = (V_tot./(R0.^2.*(1-R0./N_val).*sigma_inv)).^2;

cprob = (-2.*(Csq + Csq.^2) + lmbda.*sqrt(2.*pi).*(power(1 + Csq,2).*power(exp(1),power(lmbda,2)./2.).*erfc(lmbda./sqrt(2)) - ...
        sqrt(1 + Csq).*power(exp(1),power(lmbda,2)./(2 + 2.*Csq)).*erfc(lmbda./(sqrt(2).*sqrt(1 + Csq)))))./ ...
   (power(1 + Csq,2).*(-2 + power(exp(1),power(lmbda,2)./2.).*lmbda.*sqrt(2.*pi).*erfc(lmbda./sqrt(2))));

ind1 = 3;
ind2 = 10;

R0 = mean(R0_actual(1:ind_R0(ind1,2)-1,ind1,2));
R = 150;
phi = astar_actual2(ind1)./tot2(ind1) / alphavals(ind1);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual2(ind1)./tot2(ind1);
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
x_vals = linspace(-0.2,0.2,1000);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
S_inv_vals1 = normpdf(x_vals,0,sigma_inv);

R0 = mean(R0_actual(1:ind_R0(ind2,2)-1,ind2,2));
R = 150;
phi = astar_actual2(ind2)./tot2(ind2) / alphavals(ind2);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual2(ind2)./tot2(ind2);
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
x_vals = linspace(-0.2,0.2,1000);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;
S_inv_vals2 = normpdf(x_vals,0,sigma_inv);


figure()
h1 = histogram(fitnesses_ko(1:ind_ko(ind1,2)-1,ind1,2),'Normalization','pdf');
figure()
h2 = histogram(fitnesses_ko(1:ind_ko(ind2,2)-1,ind2,2),'Normalization','pdf');

figure()
plot((h1.BinEdges(1:end-1)+h1.BinEdges(2:end))/2,h1.Values,'k.','MarkerSize',10)
hold on
plot((h2.BinEdges(1:end-1)+h2.BinEdges(2:end))/2,h2.Values,'r.','MarkerSize',10)

%plot((h1B.BinEdges(1:end-1)+h1B.BinEdges(2:end))/2,h1B.Values,'kx','MarkerSize',2)
%hold on
%plot((h2B.BinEdges(1:end-1)+h2B.BinEdges(2:end))/2,h2B.Values,'rx','MarkerSize',2)

save('Processed_Results/FigS8C_data.mat','h1','h2','x_vals','S_inv_vals1','S_inv_vals2','-v7.3')
save('Processed_Results/FigS8D_data.mat','astar_actual1','astar_actual2','astar_actual3','tot1','tot2','tot3','coex1','coex2','coex3','alpha_star_vals','cprob','-v7.3')

close all


%% FIGURE S9A: P(Coex) for parent using R_P resources
clear all

coex2 = zeros(14,1);
tot2 = zeros(14,1);

Rp_vals = 10:10:140;

astar_actual1 = zeros(9,1);
astar_actual2 = zeros(14,1);

R0=40;
for i = 1:14
    fname = "Raw_Results/First_Step_Sims/R0_variation/sstep_astar_0.8_phi_0.1_R0_40_R_200_Rp_" + string(Rp_vals(i))+ ".mat";
    load(fname)
    for k = 1:size(deltas,1)
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/200;
        coex2(i) = coex2(i) + (d(1) > -1e-3*std(h(1-h>0))/R0)*weights(k)/mean(weights);
        tot2(i) = tot2(i) + 1;
    end
end

R0 = 40;
N_val = 200;

phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star_vals = 0.8;
alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val)*(2*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha_vals.*I_func(lmbda)*N_val/(R0*(1-R0/N_val));

eps_sq = V_tot^2 - q*R0*(1-R0/N_val)/N_val;
Rp_vals = 3:197;

V_P = sqrt(eps_sq + R0^2.*(1-Rp_vals/N_val)*q./(N_val*Rp_vals));

C = V_P./(R0*(1-Rp_vals/N_val)*sqrt(q/N_val));
L = lmbda*V_tot./V_P;

cprob = (-2.*power(C,2).*sqrt(1 + power(C,2)).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.) + ...
     L.*sqrt(2.*pi).*(power(1 + power(C,2),1.5).*power(exp(1),(1 + power(C,2)).*power(L,2)).*erfc(L./sqrt(2)) - ...
        power(exp(1),((2 + 3.*power(C,2) + 2.*power(C,4)).*power(L,2))./(2..*(1 + power(C,2)))).*erfc(L./(sqrt(2).*sqrt(1 + power(C,2))))))./...
   (power(1 + power(C,2),1.5).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.).*(-2 + power(exp(1),power(L,2)./2.).*L.*sqrt(2.*pi).*erfc(L./sqrt(2))));

save('Processed_Results/FigS9A_data','coex2','tot2','Rp_vals','cprob');

%% FIGURE S9B, S9C, S9D: Uniform distribution of R0

clear all

coex1 = zeros(10,1);
tot1 = zeros(10,1);
astar_actual1 = zeros(10,1);

coex2 = zeros(10,1);
tot2 = zeros(10,1);
astar_actual2 = zeros(10,1);

coex3 = zeros(10,1);
tot3 = zeros(10,1);
astar_actual3 = zeros(10,1);

alphavals = [0.3 0.6 1.2 2.4 4.8 9.6 19.2 38.4 76.8 153.6];

fitnesses_ko = zeros(200^2*1000,10,3);
R0_actual = zeros(200*1000,10,3);
RP_actual = zeros(200*1000,10,3);
coex_flag = zeros(200*1000,10,3);

ind_ko = ones(10,3);
ind_R0 = ones(10,3);
ind_RP = ones(10,3);


R0=(20+130)/2;
R = 150;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/R0_variation/sstep_wide_R0_variation_alpha_" + string(alphavals(i)) + "_eps_0.1_R0_20_to_130_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual1(i) = astar_actual1(i) + (length(d)-1)/R;
        coex1(i) = coex1(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot1(i) = tot1(i) + 1;

        surv_R0s = R0s{k,2};
        surv_R0s = surv_R0s(1:end-1);
        R0_actual(ind_R0(i,1):ind_R0(i,1)+length(surv_R0s)-1,i,1) = surv_R0s;
        ind_R0(i,1) = ind_R0(i,1) + length(surv_R0s);

        RP_actual(ind_RP(i,1),i,1) = surv_R0s(parent_inds(k));
        coex_flag(ind_RP(i,1),i,1) = (d(parent_inds(k)) > -1e-3*std(h)/R0);
        ind_RP(i,1) = ind_RP(i,1) + 1;

        if k <= length(organisms)
            enzymes = logical(organisms{k,1}.enzymesInSpecies);
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;

            for ind = 1:size(enzymes,1)
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,1):ind_ko(i,1)+length(inv_fitness)-1,i,1) = inv_fitness*exp(budgets(ind));
                ind_ko(i,1) = ind_ko(i,1)+length(inv_fitness);
            end
        end
    end
end


R0=(20+130)/2;
R = 150;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/R0_variation/sstep_wide_R0_variation_alpha_" + string(alphavals(i)) + "_eps_1_R0_20_to_130_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual2(i) = astar_actual2(i) + (length(d)-1)/R;
        coex2(i) = coex2(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot2(i) = tot2(i) + 1;

        surv_R0s = R0s{k,2};
        surv_R0s = surv_R0s(1:end-1);
        R0_actual(ind_R0(i,2):ind_R0(i,2)+length(surv_R0s)-1,i,2) = surv_R0s;
        ind_R0(i,2) = ind_R0(i,2) + length(surv_R0s);

        RP_actual(ind_RP(i,2),i,2) = surv_R0s(parent_inds(k));
        coex_flag(ind_RP(i,2),i,2) = (d(parent_inds(k)) > -1e-3*std(h)/R0);
        ind_RP(i,2) = ind_RP(i,2) + 1;

        if k <= length(organisms)
            enzymes = logical(organisms{k,1}.enzymesInSpecies);
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;

            for ind = 1:size(enzymes,1)
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,2):ind_ko(i,2)+length(inv_fitness)-1,i,2) = inv_fitness*exp(budgets(ind));
                ind_ko(i,2) = ind_ko(i,2)+length(inv_fitness);
            end
        end
    end
end


R0=(20+130)/2;
R = 150;
for i = 1:10
    fname = "Raw_Results/First_Step_Sims/R0_variation/sstep_wide_R0_variation_alpha_" + string(alphavals(i)) + "_eps_5_R0_20_to_130_R_150.mat";
    load(fname)
    for k = 1:length(harvests)
        h = harvests{k,1};
        h = length(h)*h/sum(h);
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        astar_actual3(i) = astar_actual3(i) + (length(d)-1)/R;
        coex3(i) = coex3(i) + (d(parent_inds(k)) > -1e-3*std(h)/R0);
        tot3(i) = tot3(i) + 1;

        surv_R0s = R0s{k,2};
        surv_R0s = surv_R0s(1:end-1);
        R0_actual(ind_R0(i,3):ind_R0(i,3)+length(surv_R0s)-1,i,3) = surv_R0s;
        ind_R0(i,3) = ind_R0(i,3) + length(surv_R0s);

        RP_actual(ind_RP(i,3),i,3) = surv_R0s(parent_inds(k));
        coex_flag(ind_RP(i,3),i,3) = (d(parent_inds(k)) > -1e-3*std(h)/R0);
        ind_RP(i,3) = ind_RP(i,3) + 1;

        if k <= length(organisms)
            enzymes = logical(organisms{k,1}.enzymesInSpecies);
            survivors = deltas{k,1} > -1e-3*std(h)/R0;
            enzymes = enzymes(survivors,:);
            delta_alphas_ko = enzymes.*(1./(sum(enzymes,2)-1)-1./sum(enzymes,2));
            alphas = enzymes.*(1./sum(enzymes,2));
            budgets = organisms{k,1}.budget;
            budgets = budgets(survivors);
        
            overall_shift = sum(alphas.*h,2);
            overall_shift = overall_shift.*exp(budgets);
            overall_shift = -mean(log(overall_shift));
            budgets = budgets + overall_shift;

            for ind = 1:size(enzymes,1)
                inv_fitness = sum(delta_alphas_ko(ind,:).*h);
                inv_fitness = inv_fitness - (delta_alphas_ko(ind,find(enzymes(ind,:),1)) + 1/(sum(enzymes(ind,:))))*h(enzymes(ind,:));
                fitnesses_ko(ind_ko(i,3):ind_ko(i,3)+length(inv_fitness)-1,i,3) = inv_fitness*exp(budgets(ind));
                ind_ko(i,3) = ind_ko(i,3)+length(inv_fitness);
            end
        end
    end
end

fitnesses_ko = fitnesses_ko(1:max(max(ind_ko))-1,:,:);

R0_actual = R0_actual(1:max(max(ind_R0))-1,:,:);
RP_actual = RP_actual(1:max(max(ind_RP))-1,:,:);
coex_flag = coex_flag(1:max(max(ind_RP))-1,:,:);

N_val = 150;
alpha_star_vals = 0.1:0.01:0.25;
R0 = mean(R0_actual(1:ind_R0(:,3)-1,:,3));
R0 = interp1(astar_actual3./tot3,R0,alpha_star_vals,'linear','extrap');
Rp_vals = mean(RP_actual(1:ind_RP(:,3)-1,:,3));
Rp_vals = interp1(astar_actual3./tot3,Rp_vals,alpha_star_vals,'linear','extrap');
phi = (astar_actual3./tot3)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
lmbda = interp1(astar_actual3./tot3,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;


alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val).*(2.*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha_vals.*I_func(lmbda).*N_val./(R0.*(1-R0./N_val));
eps_sq = V_tot.^2 - q.*R0.*(1-R0./N_val)./N_val;

V_P = sqrt(eps_sq + R0.^2.*(1-Rp_vals/N_val).*q./(N_val*Rp_vals));



C = V_P./(R0.*(1-Rp_vals./N_val).*sqrt(q/N_val));
L = lmbda.*V_tot./V_P;

cprob = (-2.*power(C,2).*sqrt(1 + power(C,2)).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.) + ...
     L.*sqrt(2.*pi).*(power(1 + power(C,2),1.5).*power(exp(1),(1 + power(C,2)).*power(L,2)).*erfc(L./sqrt(2)) - ...
        power(exp(1),((2 + 3.*power(C,2) + 2.*power(C,4)).*power(L,2))./(2..*(1 + power(C,2)))).*erfc(L./(sqrt(2).*sqrt(1 + power(C,2))))))./...
   (power(1 + power(C,2),1.5).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.).*(-2 + power(exp(1),power(L,2)./2.).*L.*sqrt(2.*pi).*erfc(L./sqrt(2))));

alpha_star_vals1 = alpha_star_vals;
cprob1 = cprob;

N_val = 150;
alpha_star_vals = 0.15:0.01:0.7;
R0 = mean(R0_actual(1:ind_R0(:,2)-1,:,2));
R0 = interp1(astar_actual2./tot2,R0,alpha_star_vals,'linear','extrap');
Rp_vals = mean(RP_actual(1:ind_RP(:,2)-1,:,2));
Rp_vals = interp1(astar_actual2./tot2,Rp_vals,alpha_star_vals,'linear','extrap');
phi = (astar_actual2./tot2)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
lmbda = interp1(astar_actual2./tot2,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;


alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val).*(2.*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha_vals.*I_func(lmbda).*N_val./(R0.*(1-R0./N_val));
eps_sq = V_tot.^2 - q.*R0.*(1-R0./N_val)./N_val;

V_P = sqrt(eps_sq + R0.^2.*(1-Rp_vals/N_val).*q./(N_val*Rp_vals));



C = V_P./(R0.*(1-Rp_vals./N_val).*sqrt(q/N_val));
L = lmbda.*V_tot./V_P;

cprob = (-2.*power(C,2).*sqrt(1 + power(C,2)).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.) + ...
     L.*sqrt(2.*pi).*(power(1 + power(C,2),1.5).*power(exp(1),(1 + power(C,2)).*power(L,2)).*erfc(L./sqrt(2)) - ...
        power(exp(1),((2 + 3.*power(C,2) + 2.*power(C,4)).*power(L,2))./(2..*(1 + power(C,2)))).*erfc(L./(sqrt(2).*sqrt(1 + power(C,2))))))./...
   (power(1 + power(C,2),1.5).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.).*(-2 + power(exp(1),power(L,2)./2.).*L.*sqrt(2.*pi).*erfc(L./sqrt(2))));

alpha_star_vals2 = alpha_star_vals;
cprob2 = cprob;


N_val = 150;
alpha_star_vals = 0.25:0.01:0.99;
R0 = mean(R0_actual(1:ind_R0(:,1)-1,:,1));
R0 = interp1(astar_actual1./tot1,R0,alpha_star_vals,'linear','extrap');
Rp_vals = mean(RP_actual(1:ind_RP(:,1)-1,:,1));
Rp_vals = interp1(astar_actual1./tot1,Rp_vals,alpha_star_vals,'linear','extrap');
phi = (astar_actual1./tot1)./alphavals';
lmbda = sqrt(2)*erfinv(1-2*phi);
lmbda = interp1(astar_actual1./tot1,lmbda,alpha_star_vals,'linear','extrap');
phi = erfc(lmbda/sqrt(2))/2;


alpha_vals = alpha_star_vals./phi;
V_tot = (1-R0/N_val).*(2.*(1 - alpha_vals.*I_func(lmbda))./(-alpha_vals .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha_vals.*I_func(lmbda).*N_val./(R0.*(1-R0./N_val));
eps_sq = V_tot.^2 - q.*R0.*(1-R0./N_val)./N_val;

V_P = sqrt(eps_sq + R0.^2.*(1-Rp_vals/N_val).*q./(N_val*Rp_vals));



C = V_P./(R0.*(1-Rp_vals./N_val).*sqrt(q/N_val));
L = lmbda.*V_tot./V_P;

cprob = (-2.*power(C,2).*sqrt(1 + power(C,2)).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.) + ...
     L.*sqrt(2.*pi).*(power(1 + power(C,2),1.5).*power(exp(1),(1 + power(C,2)).*power(L,2)).*erfc(L./sqrt(2)) - ...
        power(exp(1),((2 + 3.*power(C,2) + 2.*power(C,4)).*power(L,2))./(2..*(1 + power(C,2)))).*erfc(L./(sqrt(2).*sqrt(1 + power(C,2))))))./...
   (power(1 + power(C,2),1.5).*power(exp(1),((1 + 2.*power(C,2)).*power(L,2))./2.).*(-2 + power(exp(1),power(L,2)./2.).*L.*sqrt(2.*pi).*erfc(L./sqrt(2))));


alpha_star_vals3 = alpha_star_vals;
cprob3 = cprob;

save('Processed_Results/FigS9B_data','cprob1','cprob2','cprob3','alpha_star_vals1','alpha_star_vals2','alpha_star_vals3','astar_actual1','astar_actual2','astar_actual3','coex1','coex2','coex3','tot1','tot2','tot3');

ind1=5;

h1=histogram(R0_actual(1:ind_R0(ind1,1)-1,ind1,1),'Normalization','pdf','BinEdges',1:150);
figure()
h2=histogram(RP_actual(1:ind_RP(ind1,1)-1,ind1,1),'Normalization','pdf','BinEdges',1:150);

save('Processed_Results/FigS9C_data','h1','h2');
close all
i1 = 1;
i2 = 5;

R0 = mean(R0_actual(1:ind_R0(i1,1)-1,i1,1));
RP = mean(RP_actual(1:ind_RP(i1,1)-1,i1,1));
R = 150;
phi = astar_actual1(i1)./tot1(i1) / alphavals(i1);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual1(i1)./tot1(i1);
alpha = alpha_star/phi;

V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
x_vals = linspace(-0.05,0.05,10000);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
S_inv_vals1 = normpdf(x_vals,0,sqrt(q/R)/RP);


h1=histogram(RP_actual(1:ind_RP(i1,1)-1,i1,1),'Normalization','pdf','BinEdges',1:150);
sigma = sqrt(q/R)./(1:149)';
S_inv_vals1 = sum((1./(sigma*sqrt(2*pi))).*exp(-x_vals.^2./(2*sigma.^2)).*h1.Values');


R0 = mean(R0_actual(1:ind_R0(i2,1)-1,i2,1));
RP = mean(RP_actual(1:ind_RP(i2,1)-1,i2,1));
R = 150;
phi = astar_actual1(i2)./tot1(i2) / alphavals(i2);
lmbda = sqrt(2)*erfinv(1-2*phi);

alpha_star = astar_actual1(i2)./tot1(i2);
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
x_vals = linspace(-0.05,0.05,10000);
q = V_tot.^2.*alpha.*I_func(lmbda)*R/(R0*(1-R0/R));
S_inv_vals2 = normpdf(x_vals,0,sqrt(q/R)/RP);

h1=histogram(RP_actual(1:ind_RP(i2,1)-1,i2,1),'Normalization','pdf','BinEdges',1:150);
sigma = sqrt(q/R)./(1:149)';
S_inv_vals2 = sum((1./(sigma*sqrt(2*pi))).*exp(-x_vals.^2./(2*sigma.^2)).*h1.Values');


figure()
h1 = histogram(fitnesses_ko(1:ind_ko(i1,1)-1,i1,1),'Normalization','pdf');
figure()
h2 = histogram(fitnesses_ko(1:ind_ko(i2,1)-1,i2,1),'Normalization','pdf');

save('Processed_Results/FigS9D_data','S_inv_vals1','S_inv_vals2','x_vals','h1','h2');
close all


%% FIGURE S10: Checks for numerical consistency
clear all

astarvals = [0.5 0.8 0.99];
threshold_vals = logspace(-11,1,1000);

putative_survivors = zeros(3,1000);
tots = zeros(3,1000);

all_deltas = cell(3,2);
ext_flags = cell(3,1);
par_flags = cell(3,1);
astarvals = [0.5 0.8 0.99];


R0=40;
for i = 1:3
    fname = "Raw_Results/First_Step_Sims/knockout/sstep_astar_" + string(astarvals(i)) + "_phi_0.1_R0_40_R_200_dx_0_sigmaR2_0_muttype_-1.mat";
    load(fname)
    for k = 1:10000
        h = harvests{k,1};
        d = deltas{k,2};
        if d(end) < -1e-3*std(h(1-h>0))/R0
            continue
        end
        for t = 1:1000
            putative_survivors(i,t) = putative_survivors(i,t) + sum(deltas{k,1} > -threshold_vals(t)*std(h(1-h>0))/R0);
            tots(i,t) = tots(i,t) + 1;
        end
        survivors = deltas{k,1} > -1e-3*std(h(1-h>0))/R0;

        deltas_old = deltas{k,1};
        deltas_old = deltas_old(survivors);
        deltas_new = deltas{k,2};
        deltas_new = deltas_new(1:end-1);
        ext_flag = deltas_new < -1e-3*std(h(1-h>0))/R0;
        par_flag = (1:length(ext_flag) == parent_inds(k))';
        all_deltas{i,1} = [all_deltas{i,1}; deltas_old];
        all_deltas{i,2} = [all_deltas{i,2}; deltas_new];
        ext_flags{i} = [ext_flags{i}; ext_flag];
        par_flags{i} = [par_flags{i}; par_flag];
    end
end

save('Processed_Results/FigS10_data','threshold_vals','putative_survivors','tots','all_deltas','ext_flags','par_flags');


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

function [t_out,n,c] = log_community_dynamics(S,r,capacity,avg_supply,h_max,c_half)

    n0 = random('uniform',0,1,[S,1]);
    n0 = log((avg_supply/10)*(n0/sum(n0)));
    c0 = avg_supply*capacity/mean(capacity);

    R = length(c0);
    r0 = r;
    death_flag = false(S,1);

    [t_out,y_out] = ode45(@(t,y) f(y),[0 10^2], [n0; c0]);
    

    for i = 3:6
        y0 = y_out(end,:);
        n0 = y0(1:S);
        c0 = y0(S+1:end);

        death_flag(n0 < -35) = 1;
        n0 = n0(~death_flag);
        r0 = r(~death_flag,:);

        [tnext,ynext] = ode15s(@(t,y) f(y),[10^(i-1) 10^i],[n0 c0]);

        t_out = [t_out; tnext];

        n_next = -Inf*ones(length(tnext),S);
        n_next(:,~death_flag) = ynext(:,1:end-R);

        y_out = [y_out; [n_next ynext(:,end-R+1:end)]];
    end

    n = exp(y_out(:,1:S));
    c = y_out(:,S+1:end);

    
    function dy = f(y)
        n = y(1:end-R);
        c = y(end-R+1:end);

        dy = 0*y;
        dy(1:end-R) = sum(r0.*harvest(c,h_max,c_half)',2) - 1;
        dy(end-R+1:end) = avg_supply*capacity/mean(capacity) - harvest(c,h_max,c_half).*sum(r0.*exp(n))' - c;

        if sum(isnan(dy)) > 0
            isnan(dy)
        end
    end
end

function h = harvest(c,h_max,c_half)
    h = h_max*c./(c_half + c);
end