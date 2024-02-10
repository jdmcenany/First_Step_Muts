%% FIGURE 2B: Species sampled vs. number of survivors
clear all
close all
load('Processed_Results/Fig2B_data.mat')
figure()

p1=plot(alpha_vals_theory,astar_vals_theory(1,:),'LineWidth',1.75,'Color',[0.5 0.5 1]);
hold on
p2=plot(alpha_vals_theory,astar_vals_theory(651,:),'LineWidth',1.75,'Color',[0 0 1]);
p3=plot(alpha_vals_theory,astar_vals_theory(794,:),'LineWidth',1.75,'Color',[0 0 0.5]);
p4=plot(alpha_vals_theory,astar_vals_theory(935,:),'k','LineWidth',1.75);

errorbar(alpha_vals,astar_means(:,1),astar_stds(:,1),'.','LineWidth',1.75,'Color',[0.5 0.5 1],'MarkerSize',0.01)
errorbar(alpha_vals,astar_means(:,2),astar_stds(:,2),'.','LineWidth',1.75,'Color',[0 0 1],'MarkerSize',0.01)
errorbar(alpha_vals,astar_means(:,3),astar_stds(:,3),'.','LineWidth',1.75,'Color',[0 0 0.5],'MarkerSize',0.01)
errorbar(alpha_vals,astar_means(:,4),astar_stds(:,4),'k.','LineWidth',1.75,'MarkerSize',0.01)
xlabel('Species sampled, S/R','FontSize',25)
ylabel('Niche saturation, S*/R','FontSize',25)
legend([p1 p2 p3 p4], {'\epsilon=0.001','\epsilon=0.4','\epsilon=1.5','\epsilon=5.5'})
set(gca,'FontSize',20)
xlim([0 4])
ylim([0 1])
xticks([0 2 4])
yticks([0 0.5 1])

%% FIGURE 2C: Distribution of fitness effects
clear all
load('Processed_Results/Fig2C_data.mat')
figure()

p1=plot((edges(2:end)+edges(1:end-1))/2,h1.Values,'r.','MarkerSize',10);
hold on

p2=plot((edges(2:end)+edges(1:end-1))/2,h2.Values,'b.','MarkerSize',10);
plot(x_vals,S_inv_vals1,'k','LineWidth',1.75)
legend([p1 p2],{'Knock-Out','Knock-In'},'Location','northwest')

ylim([1 3000])
xlim([-0.01 0.01])
xlabel('Invasion fitness','FontSize',25)
ylabel('Probability density','FontSize',25)
set(gca,'YScale','log')
set(gca,'FontSize',20)

figure()

p3=plot((edges(2:end)+edges(1:end-1))/2,h3.Values,'r.','MarkerSize',10);
hold on

p4=plot((edges(2:end)+edges(1:end-1))/2,h4.Values,'b.','MarkerSize',10);
plot(x_vals,S_inv_vals2,'k','LineWidth',1.75);

legend([p3 p4],{'Knock-Out','Knock-In'},'Location','northwest')

ylim([1 3000])
xlim([-0.01 0.01])
xlabel('Invasion fitness','FontSize',25)

ylabel('Probability density','FontSize',25)
set(gca,'YScale','log')
set(gca,'FontSize',20)

%% FIGURE 2D: Sigma_inv vs. niche saturation

clear all
load('Processed_Results/Fig2D_data')
figure()

p1=plot(alpha_star,sigma_inv1*10^(3/2),'k','LineWidth',1.75);
hold on
p2=plot(alpha_star,sigma_inv2*10^(3/2),'Color',[0 0.5 0],'LineWidth',1.75);
plot(alpha_star,sigma_inv3*100^(3/2),'k','LineWidth',1.75)
plot(alpha_star,sigma_inv4*100^(3/2),'Color',[0 0.5 0])

p3=plot(astar_actual(4:end,1,1),sigma_invs(4:end,1,1)*10^(3/2),'ko','MarkerSize',6);
p4=plot(astar_actual(:,2,1),sigma_invs(:,2,1)*100^(3/2),'k.','MarkerSize',25);
plot(astar_actual(4:end,1,2),sigma_invs(4:end,1,2)*10^(3/2),'o','Color',[0 0.5 0],'MarkerSize',6)
plot(astar_actual(:,2,2),sigma_invs(:,2,2)*100^(3/2),'.','Color',[0 0.5 0],'MarkerSize',20)

set(gca,'YScale','log')
set(gca,'FontSize',20)

legend([p1 p2 p3 p4],{'S*/S = 0.1','S*/S = 0.5','R_0 = 10','R_0 = 100'})
ylim([10^-2 10^2])
xlabel('Niche saturation, S*/R','FontSize',25)
ylabel('\sigma_{inv}*R_0^{3/2}','FontSize',25)

%% FIGURE 2E: Monoculture-community fitness correlation vs. niche saturation

clear all
load('Processed_Results/Fig2E_data')
figure()

colors = [0.5 1 0.5].*(0:0.25:1)';

for i = 1:5
    plot(alpha_star,corr_theory(i,:),'Color',colors(i,:),'LineWidth',1.75)
    hold on
    plot(astar_actual(:,i),all_corrs(:,i),'k.','MarkerSize',20)
end

ylabel('Mono.-comm. correlation','FontSize',25)
xlabel('Niche saturation, S*/R','FontSize',25)
set(gca,'FontSize',20)
ylim([0 1])
xlim([0 1])
yticks([0 0.4 0.8])
xticks([0 0.5 1])
axis square

%% FIGURE 3A: P(coex) vs. niche saturation
clear all
load('Processed_Results/Fig3A_data.mat')
figure()

p1=plot(alpha_star_vals,cprob2,'Color','k','LineWidth',1.75);
hold on
p2=plot(alpha_star_vals,cprob1,'Color',[0 0.5 0],'LineWidth',1.75);
%p3=plot(alpha_star_vals,cprob1,'b','LineWidth',1.75);
%p4=plot(alpha_star_vals,cprob3,'b--','LineWidth',1.75);
p3=errorbar(astar_actual2./tot2,coex2./tot2,sqrt(coex2)./tot2,'k.','MarkerSize',10,'LineWidth',1.75);
errorbar(astar_actual1./tot1,coex1./tot1,sqrt(coex1)./tot1,'.','Color',[0 0.5 0],'MarkerSize',10,'LineWidth',1.75)


legend([p1 p2],{'S*/S = 0.1','S*/S = 0.5'},'Location','northeast')



xlabel('S*/R','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
set(gca,'FontSize',20)
set(gca,'YScale','log')
ylim([0 0.25])
ylim([0.005 0.425])
yticks([0.01 0.02 0.04 0.08 0.16 0.32])

%% FIGURE 3B: P(coex) vs. R0
clear all
load('Processed_Results/Fig3B_data.mat')
figure()
plot(R0,cprob1,'Color',[0 0 0.75],'LineWidth',2)
hold on
plot(R0,cprob2,'k','LineWidth',2)
errorbar(R0_vals,coex1./tot1, sqrt(coex1)./tot1,'.','Color',[0 0 0.75],'LineWidth',2)
errorbar(R0_vals,coex2./tot2, sqrt(coex2)./tot2,'k.','LineWidth',2)

legend({'S*/R = 0.3','S*/R = 0.8'},'Location','north')
xlabel('R_0','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
ylim([0.005 0.425])
yticks([0.01 0.02 0.04 0.08 0.16 0.32])
yticklabels([])

set(gca,'FontSize',20)
set(gca,'YScale','log')

%% FIGURE 3C: P(coex) vs. Delta X
clear all
load('Processed_Results/Fig3C_data.mat')
figure()

plot(-2:0.01:2,pcoex_vals1,'Color',[0 0 0.75],'LineWidth',1.75)
hold on
plot(-2:0.01:2,pcoex_vals2,'k','LineWidth',1.75)
errorbar(-1.5:0.3:1.5,coex1./tot1, sqrt(coex1)./tot1,'.','Color',[0 0 0.75],'LineWidth',1.75)
errorbar(-1.5:0.3:1.5,coex2./tot2, sqrt(coex2)./tot2,'k.','LineWidth',1.75)

legend({'S*/R = 0.3','S*/R = 0.8'},'Location','northeast')
xlabel('\Delta X/\sigma_{inv}','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
set(gca,'FontSize',20)
xlim([-1.6 1.6])
ylim([0.004 0.3])
set(gca,'YScale','log')
ylim([0.005 0.425])
yticks([0.01 0.02 0.04 0.08 0.16 0.32])
yticklabels([])

%% FIGURE 3C INSET: P(Coex) vs. sinv
clear all
load('Processed_Results/Fig3C_inset_data.mat')
figure()

plot(S_inv_vals,coex_theory1,'Color',[0 0 0.75],'LineWidth',1.75)
hold on
plot(S_inv_vals,coex_theory2,'k','LineWidth',1.75)
errorbar(x_vals1, cprobs(1:5,1,1),cprobs(1:5,1,2),'.','Color',[0 0 0.75],'LineWidth',1.75)
errorbar(x_vals2, cprobs(1:5,2,1),cprobs(1:5,2,2),'k.','LineWidth',1.75)


legend({'S*/R = 0.3','S*/R = 0.8'},'Location','southwest')

set(gca,'YScale','log')
xlabel('s_{inv}/\sigma_{inv}','FontSize',25)
ylabel('P(Coexistence|sinv)','FontSize',25)
ylim([10^-1 1])
xlim([0 0.5])
set(gca,'FontSize',20)

%% FIGURE 4B (with inset): Number of extinctions vs. niche saturation
clear all
load('Processed_Results/Fig4B_data.mat')
figure()
astarvals = squeeze(mean(num_species_all,2))/200;

p1=plot(astarvals(:,1), p_extinct(:,1),'k.','MarkerSize',20);
hold on
p2=plot(astarvals(:,2), p_extinct(:,2),'ko','MarkerSize',6);
legend([p2 p1],{'R0=16','R0=40'},'Location','northwest')
xlabel('S*/R','FontSize',25)
ylabel('(No. extinctions)/S*','FontSize',25)
set(gca,'FontSize',20)

figure()
errorbar(extinction_x,h1_vals/sum(h1_vals),sqrt(h1_vals)/sum(h1_vals),'k','LineWidth',1.75)
hold on
plot(extinction_x,poisspdf(extinction_x,-log(h1_vals(1)/sum(h1_vals))),'.','Color',[0.7 0.7 0.7],'MarkerSize',20)
plot(extinction_x,poisspdf(extinction_x,-log(h1_vals(1)/sum(h1_vals))),'Color',[0.7 0.7 0.7],'LineWidth',1.25)
legend({'Simulation','Poisson Fit'},'Location','northeast')
xlabel('No. Extinctions','FontSize',25)
ylabel('Probability','FontSize',25)
set(gca,'FontSize',20)
xlim([0 17])

%% FIGURE 4C: Number of resources shared between mutant and extinct strains

clear all
load('Processed_Results/Fig4C_data.mat')
figure()

ccdf1 = 1 - cumsum(h1_vals)/sum(h1_vals);
ccdf2 = 1 - cumsum(h2_vals)/sum(h2_vals);
ccdf3 = 1 - cumsum(h3_vals)/sum(h3_vals);
ccdf4 = 1 - cumsum(h4_vals)/sum(h4_vals);

p1=plot(x_bins, ccdf2,'ko','LineWidth',1.75,'MarkerSize',10);
hold on
p2=plot(x_bins, ccdf1,'k.','LineWidth',1.75,'MarkerSize',20);
p3 = plot(x_bins, ccdf4,'Color',[0.7 0.7 0.7], 'LineWidth',1.25);
plot(x_bins, ccdf3,'r','Color',[0.7 0.7 0.7], 'LineWidth',1.25);
legend({'R0=16','R0=40','Background'},'Location','northeast')

xlabel('r','FontSize',25)
ylabel('P(Shared Resources > r)','FontSize',25)
xlim([0 21])
set(gca,'FontSize',20)

%% FIGURE 4D: P(extinction) vs. relative abundance

clear all
load('Processed_Results/Fig4D_data.mat')
figure()
plot(edges(1:end-1),ext_probs,'k.')
set(gca,'XScale','log')

xlim([10^-5 10^-2.5])
xlabel('Rel. Abun.','FontSize',25)
ylabel('P(Extinction)','FontSize',25)
set(gca,'FontSize',20)

%% FIGURE 4E: Fold change in use of resource targeted by mutation in extinct species

clear all
figure()
load('Processed_Results/Fig4E_data.mat')

fold_change = zeros(10,2);
error = zeros(10,2);
for i = 1:10
    for j = 1:2
        num = totals(4,i,j)/(totals(3,i,j)+totals(4,i,j));
        den = (totals(4,i,j)+totals(2,i,j))/(totals(1,i,j)+totals(2,i,j)+totals(3,i,j)+totals(4,i,j));
        fold_change(i,j) = num/den;

        num2 = sqrt(totals(4,i,j))/(totals(3,i,j)+totals(4,i,j));
        den2 = sqrt(totals(4,i,j)+totals(2,i,j))/(totals(1,i,j)+totals(2,i,j)+totals(3,i,j)+totals(4,i,j));
        error(i,j) = fold_change(i,j)*sqrt(num2^2/num^2 + den2^2/den^2);
    end
end

p1=errorbar(R0_vals,fold_change(:,1),error(:,1),'k','LineWidth',1.75);
hold on
p2=errorbar(R0_vals,fold_change(:,2),error(:,2),'k','LineWidth',1.75);
plot([R0_vals(1) R0_vals(end)],[1 1],'k--','LineWidth',1.25)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('R_0','FontSize',25)
ylabel('Fold Change','FontSize',25)
set(gca,'FontSize',20)
ylim([10^-1.2 10^1.2])

%% FIGURE 5A, left: Example multi-step mutation trajectories
clear all
figure()

plot_ind = 133; % Highly diversifying
plot_ind = 2099; % Highly diversifying
plot_ind = 6030; % Highly diversifying


plot_ind = 4; % "Typical" surviving species
plot_ind = 488; % "Typical" extinct species
plot_ind = 86; % "Typical" extinct species
plot_ind = [760 409 168];


species0 = 1:8000;
load("Raw_Results/Multi_Step_Sims/mstep_astar_0.9_phi_0.1_R0_20_R_100_dx_0_3.mat")
last_ind = find(species_info(:,1),1);

for i = 1:(length(all_abundances)-1)
    extinct_strains = setdiff(species_present{i},species_present{i+1});

    species1 = species_present{i};
    species2 = species_present{i+1};
    abuns1 = all_abundances{i};
    abuns2 = all_abundances{i+1};
    mut_ind = find(species1 == species_info(species2(end),2));
    for j = 1:length(species1)
        original_index = find_init_species(species1(j), species_info(:,2),last_ind-1);
        ind = find(species0 == species1(j));
        if j == length(species1)
            ind = length(abuns1);
        end
        if ~isempty(find(extinct_strains == species1(j),1)) && j ~= mut_ind
            if ~isempty(find(original_index == plot_ind,1))
                plot(i,abuns1(ind),'r.','MarkerSize',20)
            else
                plot(i,abuns1(ind),'r.','MarkerSize',5)
            end
            hold on
        elseif j ~= mut_ind && ~isempty(find(original_index == plot_ind,1))
            if original_index == plot_ind(1)
                plot([i i+1],[abuns1(ind) abuns2(j)],'Color',[0 0.5 0],'LineWidth',1.25)
            elseif original_index == plot_ind(2)
                plot([i i+1],[abuns1(ind) abuns2(j)],'Color',[0 0 0.5],'LineWidth',1.25)
            else
                plot([i i+1],[abuns1(ind) abuns2(j)],'Color',[0.5 0 0],'LineWidth',1.25)
            end
            hold on
        end
    end
    original_index = find_init_species(species1(mut_ind), species_info(:,2),last_ind-1);
    ind = find(species0 == species1(mut_ind));
    if mut_ind == length(species1)
        ind = length(abuns1);
    end
    if ~isempty(find(extinct_strains == species1(mut_ind),1))
        if ~isempty(find(original_index == plot_ind,1))
            plot(i,abuns1(ind),'b.','MarkerSize',20)
        else
            plot(i,abuns1(ind),'b.','MarkerSize',5)
        end
        hold on
        if ~isempty(find(original_index == plot_ind,1))
            if original_index == plot_ind(1)
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0 0.5 0],'LineWidth',1.25)
            elseif original_index == plot_ind(2)
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0 0 0.5],'LineWidth',1.25)
            else
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0.5 0 0],'LineWidth',1.25)
            end
            hold on
        end
    else
        if ~isempty(find(original_index == plot_ind,1))
            plot(i,abuns1(ind),'Color',[0 0.9 0],'LineStyle','None','Marker','.','MarkerSize',20)
        else
            plot(i,abuns1(ind),'Color',[0 0.9 0],'LineStyle','None','Marker','.','MarkerSize',5)
        end        
        hold on
        if ~isempty(find(original_index == plot_ind,1))
            if original_index == plot_ind(1)
                plot([i i+1],[abuns1(ind) abuns2(mut_ind)],'Color',[0 0.5 0],'LineWidth',1.25)
                hold on
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0 0.5 0],'LineWidth',1.25)
            elseif original_index == plot_ind(2)
                plot([i i+1],[abuns1(ind) abuns2(mut_ind)],'Color',[0 0 0.5],'LineWidth',1.25)
                hold on
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0 0 0.5],'LineWidth',1.25)
            else
                plot([i i+1],[abuns1(ind) abuns2(mut_ind)],'Color',[0.5 0 0],'LineWidth',1.25)
                hold on
                plot([i i+1],[abuns1(ind) abuns2(length(species1)+1)],'Color',[0.5 0 0],'LineWidth',1.25)
            end
        end
    end

    species0 = species1;
end

set(gca,'YScale','log')

xlabel('No. mutations','FontSize',25)
ylabel('Rel. abun.','FontSize',25)
set(gca,'FontSize',20)
ylim([10^-5 10^-1])
xlim([1 length(all_abundances)])
xticks([1 51 101 151])
xticklabels({'0','50','100','150'})

%% FIGURE 5A, right: Abundances at which different events occur

clear all
load('Processed_Results/Fig5A_right_data.mat')
figure()
plot(h3.Values,edges,'Color',[0 0.85 0],'LineWidth',1.75)
hold on
plot(h1.Values,edges,'b','LineWidth',1.75)
plot(h2.Values,edges,'r','LineWidth',1.75)
set(gca,'YScale','log')
ylim([10^-5 10^-1])
yticks([])
xticks([])
xlim([0 max(h3.Values)*1.1])
xlabel('Frequency','FontSize',25)
legend({'Diversification','Mutation','Extinction'},'location','southeast')
set(gca,'FontSize',20)

%% FIGURE 5B, number of strains over time

clear all
load('Processed_Results/Fig5B_data.mat')
figure()


for i = 1:size(N_species,1)
    subplot(1,2,2)
    plot(N_strains(i,1:find(~N_strains(i,:),1)-1) - N_species(i,1:find(~N_species(i,:),1)-1),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
    hold on
    p2=plot(mean(N_strains(:,1:min_len)) - mean(N_species(:,1:min_len)),'Color',[0 0 0],'LineWidth',2);
    ylim([0 35])
    subplot(1,2,1)
    plot(N_strains(i,1:find(~N_strains(i,:),1)-1),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
    hold on
    p1=plot(mean(N_strains(:,1:min_len)),'Color',[0 0 0],'LineWidth',2);
    ylim([0 95])
end

subplot(1,2,1)
xlabel('No. mutations','FontSize',25)
ylabel('No. strains','FontSize',25)
set(gca,'FontSize',20)
xlim([0 200])

subplot(1,2,2)
xlabel('No. mutations','FontSize',25)
ylabel('No. close relatives','FontSize',25)
set(gca,'FontSize',20)
xlim([0 200])

%% FIGURE 5C, surival time of coexisting relatives

clear all
load('Processed_Results/Fig5C_data.mat')
figure()

plot(survival_prob,'k','LineWidth',1.75)
set(gca,'YScale','log')
xlabel('No. mutations','FontSize',25)
ylabel('P(Survival)')
set(gca,'FontSize',20)
xlim([1 95])

%% FIGURE 5D: DFE before and after evolution

clear all
load('Processed_Results/Fig5D_data.mat')
figure()

R0 = 20;
R = 100;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.9;
alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);
q = V_tot.^2.*alpha*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;

p1=plot((edges(2:end)+edges(1:end-1))/2/sigma_inv,h1.Values*sigma_inv,'r.','MarkerSize',10);
hold on

p2=plot((edges(2:end)+edges(1:end-1))/2/sigma_inv,h3.Values*sigma_inv,'b.','MarkerSize',10);
plot(x_vals/sigma_inv,S_inv_vals1*sigma_inv,'k','LineWidth',1.75)

ylim([0.5 5000]*sigma_inv)
xlim([-0.01 0.01]/sigma_inv)
legend([p1 p2],{'Before evolution','After evolution'})
xlabel('Invasion Fitness/sigmainv','FontSize',25)
set(gca,'YScale','log')
set(gca,'FontSize',20)
yticks([])

%% FIGURE S1: DFE for single organism
clear all
load('Processed_Results/Fig5E_data.mat')

figure()
means = zeros(19,1);
stds = zeros(19,1);
for i = 1:19
    coex_partial = coex_flag(:,5*(i-1)+1:5*i);
    coex_partial = coex_partial(:);

    means(i) = mean(coex_partial);
    stds(i) = std(coex_partial)/sqrt(length(coex_partial));
end


plot(2.5:5:92.5,means,'k.','MarkerSize',15)

hold on
plot(1,cprob1,'r.','MarkerSize',35)
xlabel('No. mutations','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
xlim([1 100])
legend({'Multi-step','Single-step'},'Location','southeast')
set(gca,'FontSize',20)
ylim([0 0.8])

%% FIGURE S2A: DFE for Tikhonov & Monasson consumer resource model
clear all
load('Processed_Results/FigS2A_data.mat')
figure()
p1=plot((edges(2:end)+edges(1:end-1))/2,h1.Values,'r.','MarkerSize',10);
hold on

p2=plot((edges(2:end)+edges(1:end-1))/2,h2.Values,'b.','MarkerSize',10);
plot(x_vals,S_inv_vals1,'k','LineWidth',1.75)
hold on
plot(-x_vals,S_inv_vals1,'k','LineWidth',1.75)
plot([0 0],[0.02 10^1],'k--','LineWidth',0.75)


ylim([0.02 10])
xlim([-0.2 0.2])
ylabel('Probability density','FontSize',25)
xlabel('Invasion fitness','FontSize',25)
legend([p1 p2],{'Knock-Out','Knock-In'})
set(gca,'YScale','log')
set(gca,'FontSize',20)

%% FIGURE S2B: P(coex) for knockins and knockouts in Tikhonov model
clear all
load('Processed_Results/FigS2B_data.mat')
figure()

astarvals = 0.2:0.1:0.9;
astarvals = [astarvals 0.99 0.99];
p1=plot(alpha_star_vals,pcoex_vals1,'r','LineWidth',1.75);
hold on
p2=plot(alpha_star_vals,pcoex_vals2,'b','LineWidth',1.75);
p3=plot(alpha_star_vals,pcoex_vals3,'k--','LineWidth',1.25);

astar_actual1(end-4:end-1) = [];
astar_actual2(end-4:end-1) = [];
tot1(end-4:end-1) = [];
tot2(end-4:end-1) = [];
coex1(end-4:end-1) = [];
coex2(end-4:end-1) = [];

p4=plot(astar_actual1./tot1,coex1./tot1,'kx','MarkerSize',10,'LineWidth',1.75);
p5=plot(astar_actual2./tot2,coex2./tot2,'k+','MarkerSize',10,'LineWidth',1.75);
legend([p3 p4 p5],{'Null model','Knock-Out','Knock-In'})


xlabel('S*/R','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
set(gca,'FontSize',20)
ylim([0 0.5])

%% Figure S2C: Resource availabilities for Dirichlet uptake rates

clear all
load('Processed_Results/FigS2C_data.mat')
figure()
p1=plot((edges(1:end-1)+edges(2:end))/2,h1Vals,'b.','MarkerSize',10);
hold on
p2=plot((edges(1:end-1)+edges(2:end))/2,h2Vals,'k.','MarkerSize',10);
plot(x_vals,S_inv_vals1,'k','LineWidth',1.75)
plot(x_vals,S_inv_vals2,'k','LineWidth',1.75)
set(gca,'YScale','log')
ylim([10^-2.5 10^2])
xlim([0.5 1.5])
xlabel('Resource availability','FontSize',25)
ylabel('Probability density','FontSize',25)
legend([p1 p2],{'S*/R = 0.5','S*/R = 0.8'})
set(gca,'FontSize',20)

%% Figure S2D: P(coex) for Dirichlet uptake rates

clear all
figure()
load('Processed_Results/FigS2D_data.mat')

p1=plot(alpha_star_vals,cprob,'k','LineWidth',1.75);
hold on

p3=plot(astar_actual./tot,coex./tot,'k.','MarkerSize',25);
xlabel('Niche saturation','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
set(gca,'FontSize',20)
ylim([0 0.5])

%% FIGURE S3: "Simultaneous assembly" simulations without extinction
clear all
load('Processed_Results/FigS3_data.mat')

figure()
p1=plot(alpha_star_vals,cprob,'k','LineWidth',1.75);
hold on


%p2=plot(astarvals,coex1./tot1,'k.','MarkerSize',25);
p3=errorbar(astarvals,coex1./tot1,sqrt(coex1)./tot1,'k.','LineWidth',1,'MarkerSize',25);
p4=errorbar(astarvals,coex2./tot2,sqrt(coex2)./tot2,'r.','LineWidth',0.5,'MarkerSize',20);

xlabel('S*/R','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
legend([p3 p4],{'With extinctions','Without extinctions'})
set(gca,'FontSize',20)
ylim([0 0.4])

%% FIGURE S4: P(coex) for knockout and global-effect mutations
clear all
load('Processed_Results/FigS4_left_data.mat')
figure()
subplot(1,2,1)

plot(alpha_star_vals,cprob,'k','LineWidth',1.75)
hold on
p1=plot(astar_actual1./tot1,coex1./tot1,'k.','MarkerSize',20);
p2=plot(astar_actual2./tot2,coex2./tot2,'ko','MarkerSize',7);

xlabel('Niche saturation','FontSize',25)
ylabel('P(Coexistence)','FontSize',25)
legend([p1 p2],{'Knock-out mutation','Global mutation'})
set(gca,'FontSize',20)
set(gca,'YScale','log')

ylim([0.0005 0.5])

subplot(1,2,2)
clear all
load('Processed_Results/FigS4_right_data.mat')

plot(0:0.01:1, (0:0.01:1).^2 *cprob,'k','LineWidth',1.75)
hold on
p1=plot(gammavals,coex1./tot1,'k.','MarkerSize',20);
p2=plot(gammavals,coex2./tot2,'ko','MarkerSize',7);

xlabel('Mutation effect size','FontSize',25)
legend([p1 p2],{'Knock-out mutation','Global mutation'})
set(gca,'FontSize',20)
set(gca,'YScale','log')

ylim([0.0005 0.5])
yticks([])

%% FIGURE S5: Invasion fitness vs. parent rel. abun. scatterplot
clear all
load('Processed_Results/FigS5_data.mat')
figure()

phi = 0.1;

R0 = 40;
R = 200;
phi = 0.1;
lmbda = sqrt(2)*erfinv(1-2*phi);
alpha_star = 0.8;
S = R*alpha_star/phi;

alpha = alpha_star/phi;
V_tot = (1-R0/R)*(2*(1 - alpha*I_func(lmbda))./(-alpha .*dI_func(lmbda)) - lmbda);

q = V_tot.^2.*alpha*I_func(lmbda)*R/(R0*(1-R0/R));
sigma_inv = sqrt(q/R)/R0;

slope1 = -dI_func(lmbda)*V_tot*S/(2*R0^2*(1-R0/R));
Rs = 0:0.00001:max(rel_abuns);
sinv_vals = Rs*slope1/sigma_inv;

scatter1 = scatter(s_invs(~coex_flags)/sigma_inv,rel_abuns(~coex_flags),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
hold on
scatter2 = scatter(s_invs(coex_flags)/sigma_inv,rel_abuns(coex_flags),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
scatter1.MarkerFaceAlpha = 0.05;
scatter1.MarkerEdgeAlpha = 0;
scatter2.MarkerFaceAlpha = 0.2;
scatter2.MarkerEdgeAlpha = 0;

plot(s_invs(~coex_flags)/sigma_inv,rel_abuns(~coex_flags),'Color',[0.75 0.75 0.75 0.1],'LineStyle','none','Marker','.','MarkerSize',5)
hold on
plot(s_invs(coex_flags)/sigma_inv,rel_abuns(coex_flags),'.','Color',[1 0 0 0.1],'MarkerSize',5)
plot(sinv_vals,Rs,'Color',[0.25 0.25 1],'LineWidth',2.00,'LineStyle','--')


xlabel('(Invasion fitness)/\sigma_{inv}','FontSize',25)
ylabel('Parent relative abundance','FontSize',25)
legend({'Replacement','Coexistence'})
set(gca,'FontSize',20)

axis square
ylim([0 0.05])
xlim([0 5.5])

%% FIGURE S6A, S6B: P(Coex) for non-uniform resource supply and use
clear all
figure()
load('Processed_Results/FigS6A_data');

plot(0.01:0.01:0.09,coex1./tot1,'k.','MarkerSize',20)
hold on
plot(sigmaR2_vals, cprob,'k','LineWidth',1.75)
ylim([0 0.05])
xlabel('\sigma_\kappa^2','FontSize',25)
ylabel('P(Coex)','FontSize',25)

set(gca,'FontSize',20)

clear all
figure()
load('Processed_Results/FigS6B_data');


plot(10:5:145,coex2./tot2,'k.','MarkerSize',20)
hold on
plot(Rp_vals, cprob,'k','LineWidth',1.75)
xlim([0 150])
ylim([0 0.35])
set(gca,'FontSize',20)
xlabel('R_P','FontSize',25)
ylabel('P(Coex)','FontSize',25)

%% FIGURE S7: Checks for numerical consistency
clear all
load('Processed_Results/FigS7_data');

figure()
p1=plot(threshold_vals,putative_survivors(1,:)./tots(1,:)/200,'Color',[0 0 0],'LineWidth',1.75);
hold on
p2=plot(threshold_vals,putative_survivors(2,:)./tots(2,:)/200,'Color',[0 0 0.7],'LineWidth',1.75);
p3=plot(threshold_vals,putative_survivors(3,:)./tots(3,:)/200,'Color',[0.7 0.7 1],'LineWidth',1.75);

plot([1e-3 1e-3],[0 1.5],'k--','LineWidth',0.75)
plot([1e-11 10],[0.5 0.5],'k:','LineWidth',1)
plot([1e-11 10],[0.8 0.8],'k:','LineWidth',1)
plot([1e-11 10],[0.99 0.99],'k:','LineWidth',1)


legend([p1 p2 p3],{'S*/R = 0.5','S*/R = 0.8','S*/R = 0.99'},'Location','northwest')


set(gca,'XScale','log')
xlabel('Scaled extinction threshold','FontSize',25)
ylabel('(No. putative survivors)/R','FontSize',25)
set(gca,'FontSize',20)
ylim([0 1.5])

figure()

all_deltas1 = zeros(length(all_deltas{1,1}),2);
all_deltas1(:,1) = all_deltas{1,1};
all_deltas1(:,2) = all_deltas{1,2};
all_deltas2 = zeros(length(all_deltas{2,1}),2);
all_deltas2(:,1) = all_deltas{2,1};
all_deltas2(:,2) = all_deltas{2,2};
all_deltas3 = zeros(length(all_deltas{3,1}),2);
all_deltas3(:,1) = all_deltas{3,1};
all_deltas3(:,2) = all_deltas{3,2};

subplot(1,3,1)
scatter1 = scatter(-all_deltas1(~ext_flags{1},1),-all_deltas1(~ext_flags{1},2),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
hold on
scatter2 = scatter(-all_deltas1(~~ext_flags{1},1),-all_deltas1(~~ext_flags{1},2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
scatter1.MarkerFaceAlpha = 0.01;
scatter1.MarkerEdgeAlpha = 0;
scatter2.MarkerFaceAlpha = 0.01;
scatter2.MarkerEdgeAlpha = 0;
plot([10^-16 0.1],[10^-16 0.1],'k--','LineWidth',1.75)

set(gca,'XScale','log')
set(gca,'YScale','log')

ylabel('Resource deficit post-invasion','FontSize',25)
set(gca,'FontSize',20)

xlim([10^-16 10^-5])
ylim([10^-16 0.1])
xticks([10^-14 10^-10 10^-6])
yticks([10^-15 10^-10 10^-5])

subplot(1,3,2)
scatter1 = scatter(-all_deltas2(~ext_flags{2},1),-all_deltas2(~ext_flags{2},2),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
hold on
scatter2 = scatter(-all_deltas2(~~ext_flags{2},1),-all_deltas2(~~ext_flags{2},2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
scatter1.MarkerFaceAlpha = 0.01;
scatter1.MarkerEdgeAlpha = 0;
scatter2.MarkerFaceAlpha = 0.01;
scatter2.MarkerEdgeAlpha = 0;
plot([10^-16 0.1],[10^-16 0.1],'k--','LineWidth',1.75)

set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('Resource deficit pre-invasion','FontSize',25)
set(gca,'FontSize',20)

xlim([10^-16 10^-5])
ylim([10^-16 0.1])
yticks([])
xticks([10^-14 10^-10 10^-6])


subplot(1,3,3)
scatter1 = scatter(-all_deltas3(~ext_flags{3},1),-all_deltas3(~ext_flags{3},2),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
hold on
scatter2 = scatter(-all_deltas3(~~ext_flags{3},1),-all_deltas3(~~ext_flags{3},2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); 
scatter1.MarkerFaceAlpha = 0.01;
scatter1.MarkerEdgeAlpha = 0;
scatter2.MarkerFaceAlpha = 0.01;
scatter2.MarkerEdgeAlpha = 0;
plot([10^-16 0.1],[10^-16 0.1],'k--','LineWidth',1.75)

set(gca,'XScale','log')
set(gca,'YScale','log')

set(gca,'FontSize',20)

xlim([10^-16 10^-5])
ylim([10^-16 0.1])
yticks([])
xticks([10^-14 10^-10 10^-6])

function init_species = find_init_species(ind,parentage, original_N)
    if ind <= original_N
        init_species = ind;
    else
        init_species = find_init_species(parentage(ind),parentage,original_N);
    end
end

function I = I_func(x)
    I = (1+x.^2).*erfc(x/sqrt(2))/2 - x.*exp(-x.^2/2)/sqrt(2*pi);
end

function dI = dI_func(x)
    dI = x.*erfc(x/sqrt(2)) - exp(-x.^2/2)*sqrt(2/pi);
end