function generateFigures(override)
% Generates all figures for the preprint arXiv:1609.01270:
%       M Tikhonov and R Monasson.  A collective phase in resource competition in a highly diverse ecosystem.
% 
% If override = false (default), loads the pre-computed simualtion data from disk
% If override = true, will re-compute everything. 
%
% Note, however, that each time a simulation is run for any parameters, the
% results are stored in a file named SIMDATA_REPOSITORY_*.mat so they can
% be loaded quickly the next time. 
%
% To *really* recompute everything from scratch, remove or rename that repository file.
% !! At 500 replicates of each simulation condition, this will take a long time !!
%
% Unless the goal is to exactly reproduce the figures in the paper,
% consider reducing the number of replicates, i.e. the following parameters:
% (see the code of the function "initialize")
%     params.fig2C.seedN
%     params.fig3.phaseVreplicas
%     params.fig3.phaseSreplicas
% Consider also reducing the number of elements in the following lists:
%     params.fig2C.alphaListSim       << !! if these parameters are changed,
%     params.fig3.simulationAlphaList << rerun the code with override = true !!
% Below, these parameters are marked with the comment "% SLOW" in the code.

if nargin<1
    override = false;
end

params = initialize;
data = generateData(params, override);

generate_Fig2(params, data);
generate_Fig3(params, data);
generate_SupFig_Cost(params, data);
generate_SupFig_CorrelatedCosts(params);

end

function params = initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General workflow parameters
params.saveFigs = true;
params.idStr = 'v2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display parameters
p.W = 8.6;
p.H = 4.5;
p.w = 3;
p.offX = (p.W-p.w)/2;
p.offY = 1;
params.figSXX_Cost = p;

% figure 2A
p.W = 13;
p.H = 4.2;
p.h = 2.5;
p.offY = 0.9;
p.wA = 2.5;
p.wB = 2.5;
p.wC = 2.5;
p.offXA = 0.8;
p.offXB = p.offXA+p.wA+2;
p.offXC = p.offXB+p.wB+2;
p.labelOffXA = -0.7;
p.labelOffXB = -0.7;
p.labelOffXC = -0.8;
p.labelOffY = 0.4;
params.fig2A = p;

p.W = 10.6;
p.H = 5.3;
p.wA = 3.5;
p.hA = 3.5;
p.offYA = 1.1;
p.offXA = 1.6;

p.wB = 3.3;
p.hB1 = 1;
p.hB2 = 2;
p.offXB = p.offXA+p.wA+1.8;
p.offYB2 = p.offYA;
p.offYB1 = p.offYB2+p.hA-p.hB1;
p.labelOffXA = -1.5;
p.labelOffXB = -1.4;
p.labelOffY = 0.4;
params.fig3 = p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters
K = 50;
params.fig2A.R2ListTheor = linspace(0,2,K);
params.fig2A.alphaListTheor = linspace(1,6,K);

params.fig2C.epsList = [1e-4, 0.03, 0.1];
params.fig2C.alphaListTheor = linspace(1,7,50);
params.fig2C.alphaListSim =  linspace(2,5,10);                  % SLOW
params.fig2C.seedN = 500;                                       % SLOW


params.fig3.alphaList = logspace(0,2,201);
params.fig3.phaseVreplicas = 500;                               % SLOW
params.fig3.phaseSreplicas = 500;                               % SLOW
params.fig3.simulationAlphaList = [1 2 3 4 5 6 8 10 15 20 50];  % SLOW
params.fig3.selectAlpha = [2 10];
params.fig3.simulationReplicas = 10;

clear p
p.N = 50;
p.R2 = 1;
p.sparsity = 0.5;
rFluctuationMode = [ones(p.N/2,1); -ones(p.N/2,1)];
p.capacity = 100*ones(p.N,1).*(1+sqrt(p.R2/p.N)*rFluctuationMode);
params.invariantSimulationParams = p;

end


function data = generateData(params, override)
fname = sprintf('figureData_%s.mat', params.idStr);
if ~override && exist(fname,'file')
    data = load(fname);
    data = data.data;
else   
    % Open repository of the precomputed simulation data
    repository = loadDataRepository(params);

    % Figure 2
    % Generate what's necessary to plot a phase diagram
    R2List = params.fig2A.R2ListTheor;
    alphaList = params.fig2A.alphaListTheor;

    p = params.invariantSimulationParams;
    p.epsilon = 1e-7;
    psi = NaN(length(alphaList), length(R2List));
    for i=1:length(alphaList)
        for j=1:length(R2List)
            p.alpha = alphaList(i);
            p.R2 = R2List(j);
            try
                par = solveForParams(p);
                psi(i,j) = par.psi;
            end
        end
    end
    d.psi = psi;
    d.R2List = R2List;
    d.alphaList = alphaList;
    data.fig2A = d;


    epsList = params.fig2C.epsList;
    alphaList = params.fig2C.alphaListTheor;
    p = params.invariantSimulationParams;
    survTheor = zeros(length(epsList),length(alphaList));
    for e=1:length(epsList)
        p.epsilon = epsList(e);
        for a=1:length(alphaList)
            p.alpha = alphaList(a);
            par = solveForParams(p);
            survTheor(e,a) = par.survivorFraction;
        end
    end
    d.survTheor = survTheor;
    d.alphaListTheor = alphaList;

    alphaList = params.fig2C.alphaListSim;
    seedN = params.fig2C.seedN;
    p = params.invariantSimulationParams;
    survSim = zeros(length(epsList),length(alphaList),seedN);
    for e=1:length(epsList)
        for a=1:length(alphaList)
            for s=1:seedN
                fprintf('alpha = %f, seed %d of %d.\n',alphaList(a),s, seedN);
                [repository,h,prop] = computeOrLoad(repository, alphaList(a), s, epsList(e));
                delta = (prop.enzymesInSpecies*h(:) - prop.cost)./prop.enzCount;
                % Record the number of survivors
                deltaSrt = sort(delta, 'descend'); thresh = deltaSrt(p.N);
                surv = delta >= max(-1e-8, thresh);
                survSim(e,a,s) = sum(surv);
            end
        end
    end
    d.survSim = survSim/p.N;
    d.alphaListSim = alphaList;
    d.epsList = epsList;
    data.fig2C = d;

    
    % Figure 3A: crossover
    alphaList = params.fig3.alphaList;
    p = params.invariantSimulationParams;
    % preallocate memory; parameter values irrelevant
    p.epsilon = 1e-5;  p.alpha = 50;
    pList(length(alphaList)) = solveForParams(p);

    % theory data
    for power = 1:5
        p.epsilon = 10^(-power);
        for a=1:length(alphaList)
            p.alpha = alphaList(a);
            pList(a) = solveForParams(p);
        end
        fieldName = sprintf('pList%d',power);
        d.(fieldName) = pList;       
    end
    
    % simulation data
    p = params.invariantSimulationParams;
    % For a range of alpha, and for each, repeat for a number of tries
    % Run the simulation for a range of epsilon, recording harvests
    maxPower = 5;
    seedN = params.fig3.simulationReplicas;
    harvestTbl = zeros(maxPower,length(params.fig3.simulationAlphaList), seedN, p.N);
    for alphaInd = 1:length(params.fig3.simulationAlphaList)
        % check if data for this alpha was already computed
        p.alpha = params.fig3.simulationAlphaList(alphaInd);
        % Now for each seed, run convergence for a range of epsilon
        for seed = 1:seedN
            fprintf('alpha = %f, seed %d of %d.\n',p.alpha,seed, seedN);
            for power = 1:maxPower
                [repository, harvest] = computeOrLoad(repository, p.alpha, seed, 10^(-power));
                harvestTbl(power, alphaInd, seed,:) = harvest;
            end
        end 
    end
     
    d.alphaListTheory = alphaList; 
    d.alphaListSimulation = params.fig3.simulationAlphaList;
    d.harvestTable = harvestTbl; 
    data.fig3A = d;
    clear d

    % Figure 4B: simulations with extra copies, for panel B
    epsilon = 1e-3;
    [repository, d.infoV] = computeDetailedPhaseInfo(repository, params.fig3.selectAlpha(1), epsilon, params.fig3.phaseVreplicas);
    [repository, d.infoS, deltaTable] = computeDetailedPhaseInfo(repository, params.fig3.selectAlpha(2), epsilon, params.fig3.phaseSreplicas); 
    data.fig3B = d;
    
    % Use the S-phase data for the histogram in fig. 2B
    data.fig2B.deltaTable = deltaTable;
    data.fig2B.params = params.invariantSimulationParams;
    data.fig2B.params.epsilon = epsilon;
    data.fig2B.params.alpha = params.fig3.selectAlpha(2);
    
    saveDataRepository(repository);

    save(fname, 'data');
end
end

function [repository, info, deltas] = computeDetailedPhaseInfo(repository, alpha, epsilon, seedN)
    p = repository.invariantParams;
    p.alpha = alpha;
    p.epsilon = epsilon;
    totalSpecies = p.alpha*p.N;

    harvest = zeros(seedN, p.N);
    performanceExtEnv = NaN(seedN, totalSpecies);
    success = NaN(seedN, totalSpecies);
    cpe = NaN(seedN, totalSpecies);
    survCount = zeros(seedN, 1);
    if nargout>2
        deltas = zeros(seedN, round((p.alpha-1)*p.N));
    end
    for seed=1:seedN
        fprintf('alpha = %f, seed %d of %d.\n',p.alpha,seed, seedN);
        [repository, h,prop] = computeOrLoad(repository, p.alpha, seed, p.epsilon);
        harvest(seed,:) = h;
        
        % now compute (for ALL species, not just survivors):
        % their "success", cost per enzyme, and performance in ext. env.
        
        % cost:
        cpe(seed,:) = prop.cost./prop.enzCount;
        
        % performance in external environment:
        extEnv = prop.capacity;
        normalizedEnzTable = prop.enzymesInSpecies./repmat(prop.enzCount(:),[1, p.N]);
        performanceExtEnv(seed,:) = (normalizedEnzTable*extEnv)' - cpe(seed,:);
        
        % success:
%        delta = (prop.enzymesInSpecies*h(:) - prop.cost)./prop.enzCount;
        delta = (prop.enzymesInSpecies*h(:) - prop.cost); % Do not renormalize the strategy vectors 
        
        % Record the number of survivors
        [deltaSrt,order] = sort(delta, 'descend'); thresh = deltaSrt(p.N);
        surv = delta >= max(-1e-8, thresh);
        survCount(seed) = sum(surv);
        % save the deltas if we need to
        if nargout>2
            deltas(seed,:) = deltaSrt(p.N+1:end);
        end
        
        % Find abundance of survivors at equilibrium:
        abdSurv = (prop.capacity'./h)*pinv(double(prop.enzymesInSpecies(surv,:)));
        abdSurv = max(abdSurv, 0.1); % eliminate negatives
        abd0 = zeros(size(prop.cost)); abd0(surv)=abdSurv;
        % to verify, run convergence just in case 
        % (but in fact abd1 and abd0 should be almost identical)
        abd1 = equilibrateCommunity(abd0,prop,1e-9);        
        
        % success (at true equilibrium)
        % for survivors, success = their abundance (positive) 
        % for others, their delta (negative)
        s = abd1;
        s(abd1==0) = delta(abd1==0);
        success(seed,:) = s;
    end
    info.harvest = harvest; % the V-phase
    info.performanceExtEnv = performanceExtEnv;
    info.success = success;
    info.cpe = cpe;
    info.survCount = survCount;

    info.capacity = p.capacity;
    info.totalSpeciesCount = p.alpha * p.N;
end

% a respoitory stores the resutls of calculations with some parameters
% fixed, while others are changing. When creating a new repository, the 
% invariant parameters are copied directly from
% params.invariantSimulationParams. When loading an existing one,
% consistency is checked. In all simulations, alpha, seed and espilon are
% supplied, all other parameters are loaded dircetly from invariantParams.
function repository = loadDataRepository(params)
fname = sprintf('SIMDATA_REPOSITORY_%s.mat', params.idStr);
if exist(fname,'file')
    repository = load(fname);
    repository = repository.repository;
    % verify consistency of global repository parameters
    assert(isequal(repository.invariantParams, params.invariantSimulationParams),...
        'Parameters values stored in the repository are inconsistent with this simulation run; repository must be recomputed.');
else
    repository.fname = fname;
    repository.invariantParams = params.invariantSimulationParams;
    % repository format: 
    % headers array that contains changeable parameters in the order 
    %   alpha, seed, epsilon
    % harvest table that contains the computed harvest values
    repository.entries = 0;
    % preallocate memory for 1000 entries
    repository.header = NaN(1000,3);
    repository.harvest = NaN(1000,repository.invariantParams.N);
end
end

function saveDataRepository(repository)
save(repository.fname,'repository');
end

% yes one could do something much smarter with hashing, but this is a
% one-time operation and so I don't care that much about efficiency here
function harvest = retrieveFromRepository(repository, header)
compare = all(repository.header(1:repository.entries, :) == repmat(header,[repository.entries, 1]),2);
index = find(compare,1);
harvest = repository.harvest(index,:); % will be empty if no matches were found
end

function repository = saveToRepository(repository, header, harvest)
if repository.entries==size(repository.header,1)
    % repository full. Allocate a new chunk.
    repository.header(end+1000,1)=NaN;
    repository.harvest(end+1000,1)=NaN;
end
repository.entries = repository.entries+1;
repository.header(repository.entries,:) = header;
repository.harvest(repository.entries,:) = harvest;
if mod(repository.entries,10) == 0
    % save to disk with every 10 new entries
    saveDataRepository(repository);
end
end

function [repository, harvest, prop] = computeOrLoad(repository, alpha, seed, epsilon)
    % load pre-computed result if it exists
    harvest = retrieveFromRepository(repository, [alpha, seed, epsilon]);
    % if not, need to compute it
    if isempty(harvest)
        prop = seed2prop(repository, alpha, seed, epsilon);
        harvest = findEquilibrium(prop);
        repository = saveToRepository(repository, [alpha, seed, epsilon], harvest);        
    else
        % harvests are loaded from file, so we have that part - but not
        % prop. If requested to return prop as an output, need to generate it
        if nargout>2            
            prop = seed2prop(repository, alpha, seed, epsilon);
        end
    end
end

function prop = seed2prop(repository, alpha,seed,epsilon)
    p = repository.invariantParams;
    p.alpha = alpha;
    p.seed = round(1/epsilon+seed);
    p.epsilon = epsilon;
    prop = generateOrganisms(p);
    prop.capacity = p.capacity;
end


function generate_SupFig_Cost(params, data)
p = params.figSXX_Cost;
W = p.W;
H = p.H;
w = p.w;

clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
    'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
    'Units','Centimeters','Position',[14 7 W H]); 

%%% PANEL A %%%
axes('Units','Centimeters',...
    'Position',[p.offX, p.offY, w, w]);

x=linspace(0,1,20);
rng(0)

chi0 = 1;
dMean = sin(pi*x);
dChi = random('norm', dMean, sin(pi*x)/2);
chi = chi0 * (1+dChi);
%clf;
hold on;
plot(x,chi0*(1+dMean),'r-')
plot(x,chi,'k.','MarkerSize',7)
sel = chi<chi0+1e-6;
plot(x(sel),chi(sel),'k.','MarkerSize',20)
plot([0 1],chi0*[1 1], 'k--')
axis([0 1 0 2.5])
text(0.6,1,'\it\chi = \chi_0','VerticalAlignment','top','FontSize',12)
text(0.5,-0.2,'\itx','VerticalAlignment','top','HorizontalAlignment','center','FontSize',16)
%annotation('textarrow',0.1+[0.5627 0.4984],[0.7619 0.7178],'String','{\it\chi}/{\it\chi_0} = 1 + sin {\it\pix}',...
%    'HeadLength',8,'HeadWidth',8,'FontSize',12)
ylabel('Cost')
set(gca,'YTick',0:2,'YTickLabel',{'0','\chi_0','2\chi_0'},'XTick',0:1);
adjustSizes(gca,1,14);
%%
if params.saveFigs
    fname = sprintf('FigureSXX_CostModel_%s', params.idStr);
    saveas(gcf,[fname '.pdf'])
    saveas(gcf,[fname '.png'])
    saveas(gcf,[fname '.fig'])
end
end

function generate_Fig2(params, data)
p = params.fig2A;
W = p.W;
H = p.H;
%%
clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
    'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
    'Units','Centimeters','Position',[14 7 W H]); 

%%% PANEL A %%%
d = data.fig2A;
axes('Units','Centimeters',...
    'Position',[p.offXA, p.offY, p.wA, p.h]);
logPsi = log(d.psi);
logPsi(logPsi<-8)=-Inf;
im = imagesc(d.R2List([1,end]),d.alphaList([1,end]),logPsi);
axis([0 2 1 5.2])
a=axis;
set(gca,'YTick',[2 4],'XTick',[0 1 2],'XTickLabel',{'0','','2'},'TickLength',[0.04 0.1]);
text(mean(a([1,2])),a(3)-0.4,'$\overline{\delta R^2}$','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalALignment','top','FOntSize',12);
text(-0.35,mean(a([3,4])),'$\alpha$','Interpreter','latex',...
    'HorizontalAlignment','right','VerticalALignment','middle','FOntSize',12);
axis xy
axis square
title('$\log \psi$','Interpreter','latex','FontSize',12)
colormap(redblue)

% Overlay the parameteric prediction for the critical line
mus=0:0.001:2;
alphaLine = 1./HH(mus);
sparsity=0.5;
dR2Line = mus.^2./(1-II(mus)./HH(mus))*(1-sparsity)/sparsity;
hold on;
% plot(dR2Line,alphaLine,'c-','LineWidth',2)
% plot(dR2Line,alphaLine,'m--','LineWidth',2)
%plot(dR2Line,alphaLine,'g-','LineWidth',2)
plot(dR2Line,alphaLine,'k:','LineWidth',2)

[~, muCrit] = min(abs(dR2Line-1));
alphaCritical = alphaLine(muCrit);

[figx, figy] = dsxy2figxy(gca, [1 1], [2 5]);
annotationColor = 'k'; %[0 0.6 0]; %0.3*[1 1 1];
arrowStyle = {'HeadLength',8,'HeadWidth',8,'HeadStyle','vback1','Color',annotationColor,'LineWidth',2.5};
annotation('arrow',figx,figy,arrowStyle{:}); 
hold on
plot(1, alphaCritical,'.','Color',annotationColor,'MarkerSize',25); 

cb = colorbar('Location','manual','Units','centimeters','Position',[p.offXA+p.wA+0.3, p.offY, 0.3, p.h]);
cb.Limits(1) = -8;
cb.TickLabels{1}='-\infty';

% wash out the image 
im.AlphaData = 0.5;
% wash out the colorbar accordingly
annotation('textbox',...
    'units','centimeters','Position',cb.Position,...
    'FitBoxToText','off',...
    'FaceAlpha',1-im.AlphaData,...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

text(p.labelOffXA,p.h+p.labelOffY,'(a)', 'Units', 'Centimeters', 'FontSize', 16);
text(0.2,5,'S','FontSize',18,'Color','k','VerticalAlignment','top','FontWeight','bold');
text(1.4,1.9,'V','FontSize',18,'Color','k','FontWeight','bold');

% 
% 
% 
% lim = 1.2;
% hold on
% plot([0 1], [1 1], 'k:', 'LineWidth',1);
% plot([1 1], [0 1], 'k:', 'LineWidth',1);
% x = [0.3 0.8 0.9];
% y = [0.9 0.75 0.5];
% patch([0 0 x(1) x(2) x(3) x(3)],[0 y(1) y(1) y(2) y(3) 0],0.6*[1 1 1]);
% red = [0.8 0 0];
% %plot(x(2),y(2),'k.','MarkerSize',10);
% plot(1,1,'.','MarkerSize',10,'Color',red);
% text(0.42,0.42,'\Omega','FontSize',24,'HorizontalAlignment','center','VerticalALignment','middle');
% axis([0 lim 0 lim]);
% set(gca,'XTick',[], 'YTick', []);
% %set(gca,'XTick',[0 1], 'YTick', [0 1], 'XTickLabel',{'0','\chi_0'}, 'YTickLabel',{'0','\chi_0'})
% adjustSizes(gca,1,12);
% [figx, figy] = dsxy2figxy(gca, [1 x(2)], [1 y(2)]); 
% annotation('arrow',figx,figy,'HeadLength',8,'HeadWidth',8,'HeadStyle','vback1','Color',red); 
% text(1.03,1,'$\vec g$','Color',red,'HorizontalAlignment','left','VerticalALignment','top','Interpreter','Latex','FontSize',18);
% % text(0.05,1.1,'\it\chi\rm_0', 'FontSize', 12,'HorizontalAlignment','left','VerticalALignment','baseline');
% % text(1.02, 0,'\it\chi\rm_0', 'FontSize', 12,'HorizontalAlignment','left','VerticalALignment','bottom');
% text(1.1, 1.03,'$(\chi_0,\chi_0)$', 'FontSize', 12,'HorizontalAlignment','right','VerticalALignment','bottom','Interpreter','latex');
% text(0.5,-0.25,'\ith\rm_1', 'FontSize',12);
% text(-0.3,0.5,'\ith\rm_2', 'FontSize',12);
% text(p.labelOffXA,p.h+p.labelOffY,'A', 'Units', 'Centimeters', 'FontSize', 20);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = data.fig2B;
FACTOR = 1e3; % to remove the annoying scaling factor from x axis
axes('Units','Centimeters',...
    'Position',[p.offXB, p.offY, p.wB, p.h]);
theory = solveForParams(d.params);
%ctr = theory.mu * theory.psi / (d.params.sparsity*d.params.N)*FACTOR;
ctr = theory.mu * theory.psi *FACTOR;
step = 0.04;
binctr = (-5:step:0)*ctr;
%distrib = pdf('normal',binctr, -ctr, theory.psi/(d.params.sparsity*d.params.N)*FACTOR)*step*ctr;
distrib = pdf('normal',binctr, -ctr, theory.psi*FACTOR)*step*ctr;
% extra = 5;
% binctrExtended = [binctr,(1:extra+1)*step*ctr]';
% distrib(end:(end+extra)) = (1-sum(distrib(1:end-1)))/(extra+1);
% distrib(end+1)=0;
cts = histcounts(d.deltaTable(:)*FACTOR,binctr)/(d.params.alpha*d.params.N*size(d.deltaTable,1));
%bar(binctr(1:end-1)+diff(binctr),cts,'FaceColor',0.8*[1 1 1],'EdgeColor',0.8*[1 1 1],'LineWidth',1);
hold on
deltaPeakHeight = 2.2e-2;
patch([-2 -2 2 2]*step*ctr,[0 deltaPeakHeight*[1 1] 0],'k','EdgeColor','none')
plot(binctr, distrib,'k-','LineWidth',2)
plot(binctr(1:end-1)+diff(binctr),cts,'r.');
x0 = -ctr; y0 = distrib(binctr==x0);
plot([x0,x0],[0,y0],'k--','LineWidth',1)
text(x0,-3e-4,'$pm$','Interpreter','LaTeX','HorizontalAlignment','center','VerticalAlignment','top','FontSize',12);
axis([-6.25e-3*FACTOR,1e-3*FACTOR,0,2.5e-2]);
set(gca,'YTick',[],'YAxisLocation','origin','box','off','XTick',[-5 0],'XTickLabel',{'$-5\epsilon$','0'},'TickLabelInterpreter','LaTeX');
[figx, figy] = dsxy2figxy(gca, [5 4]*step*ctr, deltaPeakHeight*[1 1]); 
text(16*step*ctr,deltaPeakHeight,'$E(\lambda)$','Interpreter','LaTeX','FontSize',10,'VerticalAlignment','middle')
annotation('arrow',figx,figy,'HeadLength',8,'HeadWidth',8,'HeadStyle','vback1'); 

selX = [78,124];
[figx, figy] = dsxy2figxy(gca, binctr(selX), distrib(selX)); 
annotation('doublearrow',figx,figy,'Head1Length',8,'Head1Width',8,'Head1Style','vback1','Head2Length',8,'Head2Width',8,'Head2Style','vback1'); 
text(mean(binctr(selX)),0.9*mean(distrib(selX)),'$\psi$','Interpreter','LaTeX','FontSize',12,'VerticalAlignment','top','HorizontalAlignment','center','Background','w')

text(1.5,-0.1e-2,'$\Delta$','Interpreter','LaTeX','FontSize',12,'HorizontalAlignment','left')
adjustSizes(gca,1,10)
title('$p(\Delta)$','Interpreter','LaTeX','FontSize',12)
text(p.labelOffXB,p.h+p.labelOffY,'(b)', 'Units', 'Centimeters', 'FontSize', 16);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = data.fig2C;
axes('Units','Centimeters',...
    'Position',[p.offXC, p.offY, p.wC, p.h]);

hold on
patch(alphaCritical*[0 0 1 1], 1.1*[0 1 1 0],'r','FaceAlpha',0.1,'EdgeColor','none');
patch(alphaCritical*[1 1 2 2], 1.1*[0 1 1 0],'b','FaceAlpha',0.1,'EdgeColor','none');
plot(alphaCritical*[1 1],[0 1.1],'k:','LineWidth',2)
% Use the simulation data to infer error bars for the plot in fg2C
% Procedure: for each alpha used in simulation, compute the std of
% log(q). Then interpolate these values for all other alpha
for i=1:3
    mn = squeeze(mean(d.survSim(i,:,:),3));
    sd = squeeze(std(d.survSim(i,:,:),[],3));
    ste = sd; 
    ste = 1*sd./sqrt(size(d.survSim,3));
    lower = mn-ste;
    upper = mn+ste;
    xs = params.fig2C.alphaListSim;
    patch([xs, fliplr(xs)],[upper, fliplr(lower)],'k','FaceAlpha',0.2,'EdgeColor','none');
    %patch([xs; xs(end:-1:1)],[upper; lower(end:-1:1)],'k','FaceAlpha',0.2,'EdgeColor','none');
    plot(xs, mn, 'r-');
end
plot(d.alphaListTheor,d.survTheor(1,:),'k-','LineWidth',1)
plot(d.alphaListTheor,d.survTheor(2,:),'k-','LineWidth',1)
plot(d.alphaListTheor,d.survTheor(3,:),'k-','LineWidth',1)

axis([2 5 0.74 1.003])
[figx, figy] = dsxy2figxy(gca, [2 5], 0.74*[1 1]);
annotation('arrow',figx,figy,arrowStyle{:}); 
plot(alphaCritical, 0.74, '.','Color',annotationColor,'MarkerSize',25)

% mn1 = squeeze(mean(d.survSim(1,:,:),3));
% mn2 = squeeze(mean(d.survSim(2,:,:),3));
% mn3 = squeeze(mean(d.survSim(3,:,:),3));
% plot(d.alphaListSim, mn1, 'r-');
% plot(d.alphaListSim, mn2, 'r-');
% plot(d.alphaListSim, mn3, 'r-');
% plot(alphaSimTable(:),squeeze(d.survSim(1,:))','r.')
% plot(alphaSimTable(:),squeeze(d.survSim(2,:))','g.')
% plot(alphaSimTable(:),squeeze(d.survSim(3,:))','b.')
text(4.9,1,'1e-4','VerticalAlignment','top','HorizontalAlignment','right');
text(4.9,0.925,'0.03','VerticalAlignment','top','HorizontalAlignment','right');
text(4.9,0.84,'0.1','VerticalAlignment','top','HorizontalAlignment','right');
%plot(alphaCritical*[1 1],[0,1],'r-')
title('Survivors','FontSize',12,'FontWeight','normal')
set(gca,'XTick',[2:5],'YTick',0.8:0.1:1,'YTickLabel',{'$0.8\,N$', '$0.9\,N$', '$N$'},'TickLabelInterpreter','LaTeX')
text(3.8,0.74,'S','FontSize',12,'Color','b','FontWeight','bold','VerticalAlignment','bottom');
text(3.24,0.74,'V','FontSize',12,'Color','r','FontWeight','bold','VerticalAlignment','bottom');
text(3.5,0.68,'$\alpha$','Interpreter','LaTeX','FontSize',12,'HorizontalAlignment','center');
text(p.labelOffXC-0.4,p.h+p.labelOffY,'(c)', 'Units', 'Centimeters', 'FontSize', 16);




%%
if params.saveFigs
    fname = sprintf('Figure2_%s', params.idStr);
    saveas(gcf,[fname '.pdf'])
    saveas(gcf,[fname '.fig'])
    saveas(gcf,[fname '.png'])
end
end


function generate_Fig3(params, data)
p = params.fig3;
W = p.W;
H = p.H;
%%
maxPower = 5;
clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
    'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
    'Units','Centimeters','Position',[14 7 W H]); 

%%% PANEL A %%%
d = data.fig3A;

axes('Units','Centimeters',...
    'Position',[p.offXA, p.offYA, p.wA, p.hA]);

loglog(d.alphaListTheory, [d.pList1.psi],'k-',d.alphaListTheory,[d.pList2.psi],'k-',  d.alphaListTheory,[d.pList4.psi],'k-',d.alphaListTheory,[d.pList5.psi],'k-', ...
    'LineWidth',1);
hold on;
loglog(d.alphaListTheory, [d.pList3.psi],'k-', 'LineWidth',1);
%title('\itq') axes
axis([0.8 7e2 8e-6 1])
%set(gca,'YTick',[1e-10 1e-5 1],'YTickLabel',{'10^{-10}', '10^{-5}', '1'},'XTick',[1 10 100],'YScale','log','XScale','log')
set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 0.1 1],'YTickLabel',{'10^{-6}', '', '10^{-4}', '', '10^{-2}', '', '1'},...
    'XTick',[1 10 100],'YScale','log','XScale','log','YGrid','on','YMinorGrid','off')
adjustSizes(gca,1,12);
label = {' \fontsize{12}\epsilon\fontsize{10} = 0.1','0.01','1e-3','1e-4','1e-5'};

% % When interpolating the error bars, use this vector of alphas:
% maxAlpha = max(d.alphaListSimulation);
% alphaVec = [d.alphaListTheory, maxAlpha];
% alphaVec = unique(alphaVec(alphaVec<=maxAlpha)); % this includes all alpha used for plotting the theory curve, up to and including maxAlpha
%
for i=1:maxPower
    fieldName = sprintf('pList%d',i);
    pList = d.(fieldName);
    yCoord = pList(round(end)).psi;
    text(6e2,yCoord,label{i},'HorizontalAlignment','right','VerticalAlignment','baseline','fontSize',10,'BackgroundColor','w','Margin',1);
    % Show datapoints from simualtions
    alphaTbl = repmat(p.simulationAlphaList,[size(d.harvestTable,3),1])';
    g=1-d.harvestTable(i, :, :, :);
    qTbl = squeeze(sum(g.^2,4));
    sp = pList(1).sparsity;
    epsilon = pList(1).epsilon;
    psiTbl = sqrt(qTbl(:)*sp*(1-sp)+epsilon^2);
    plot(alphaTbl(:), psiTbl(:),'.')

%     % draw the shaded area based on simulation results
%     fieldName = sprintf('simData%d',i);
%     stdInterp = interp1(d.alphaListSimulation, d.(fieldName).stdLog,alphaVec);
%     meanInterp = interp1(d.alphaListSimulation, d.(fieldName).meanLog,alphaVec);
%     lower = exp(meanInterp+stdInterp);
%     upper = exp(meanInterp-stdInterp);
%     patch([alphaVec, fliplr(alphaVec)], [upper,fliplr(lower)],0.8*[1 1 1],'EdgeColor','none');
    
end
val = interp1(d.alphaListTheory, [d.pList3.psi], p.selectAlpha);
plot(p.selectAlpha(1),val(1),'r.','MarkerSize',20);
text(p.selectAlpha(1),val(1),'V','Color','r','FontSize',16,'VerticalAlignment','top','HorizontalAlignment','right','FontWeight','bold');
plot(p.selectAlpha(2),val(2),'b.','MarkerSize',20);
text(p.selectAlpha(2)-2,val(2)*1.4,' S','Color','b','FontSize',16,'VerticalAlignment','baseline','HorizontalAlignment','left','FontWeight','bold');
%text();
text(30, 6e-7, '\alpha', 'FontSize', 16, 'VerticalAlignment','middle','HorizontalAlignment','center')
text(p.labelOffXA,p.hA+p.labelOffY,'(a)', 'FontSize', 16, 'units','centimeters');
%text(0.08,1e-3,'fluctuations of \ith_i', 'FontSize', 12, 'VerticalAlignment','middle','HorizontalAlignment','center','rotation',90);
text(0.08,1e-3,'\psi', 'FontSize', 14, 'VerticalAlignment','middle','HorizontalAlignment','center');
%
%%% PANEL B, top %%%
d = data.fig3B;
axes('Units','Centimeters',...
    'Position',[p.offXB, p.offYB1, p.wB, p.hB1]);
resourceSupply = d.infoV.capacity;
plot(resourceSupply./mean(resourceSupply),'k.-');
a=[1 50 0.7 1.3];
axis(a);
set(gca,'XTickLabel',[],'YTick',a(3:4));
title('External supply','FontSize',12,'FontWeight','normal');
text(-11,1,'\itR_i','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
box on

%%% PANEL B, bottom %%%
axes('Units','Centimeters',...
    'Position',[p.offXB, p.offYB2, p.wB, p.hB2]);
set(gca,'XTick',[1,50]);
text(25.5, 0.825, 'resource','FontSize',12,'HorizontalAlignment','center');
text(-11,1,'\ith_i','HorizontalAlignment','center','FontSize',12);

mn = mean(d.infoV.harvest);
sd = std(d.infoV.harvest);
xs = 1:length(mn);

upper = mn+sd;
lower = mn-sd;

hold on
h1 = plot(xs,mn,'r-','LineWidth',2);
patch([xs, fliplr(xs)],[upper, fliplr(lower)],'r','FaceAlpha',0.2,'EdgeColor','none');
text(26,mean(upper(1:25)),'V-phase','Color','r','FontSize',10);

mn = mean(d.infoS.harvest);
sd = std(d.infoS.harvest);
xs = 1:length(mn);
upper = mn+sd;
lower = mn-sd;
patch([xs, fliplr(xs)],[upper, fliplr(lower)],'b','FaceAlpha',0.2,'EdgeColor','none');
h1 = plot(xs,mn,'b-','LineWidth',2);
text(2,mean(lower(1:25)),'S-phase','Color','b','FontSize',10,'VerticalAlignment','top');

text(p.labelOffXB,-p.offYB2+p.offYA+p.hA+p.labelOffY,'(b)', 'Units', 'Centimeters', 'FontSize', 16);
axis([1 50 0.9 1.1])
box on

if params.saveFigs
    fname = sprintf('Figure3_%s', params.idStr);
    saveas(gcf,[fname '.pdf'])
    saveas(gcf,[fname '.png'])
    saveas(gcf,[fname '.fig'])
end

%% Sup figure
clf;
W = 17;
H = 16;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
    'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
    'Units','Centimeters','Position',[1 7 W H]); 
% Scatter plot performance rank in external environment vs abundance at equilibrium
[cVec, cVecBaseline] = evaluatePredictivePower(d.infoS);

subplot(2,2,3);
sel = d.infoS.success(1,:)>0; %<=d.infoS.survCount(1);
rkAbd = getRank(-d.infoS.success(1,:));
rkPerf = getRank(-d.infoS.performanceExtEnv(1,:));
plot(rkPerf(~sel),rkAbd(~sel),'.','MarkerSize',15,'Color',[0.6 0.6 1])
hold on;
plot(rkPerf(sel),rkAbd(sel),'.','MarkerSize',15,'Color',[0 0 1])
xlabel('Performance rank when alone');
ylabel('Actual success rank');
title(sprintf('S-phase example; c = %.2f ',cVec(1)));
axis([0 d.infoS.totalSpeciesCount 0 d.infoS.totalSpeciesCount]);
dx=8;
off = -1.7;
dh = 6;
opts = {'FontSize',18,'Units','Centimeters','BackgroundColor','w','Margin',5};
text(off,dh,'(c)',opts{:})
text(off+dx,dh,'(d)',opts{:})

subplot(2,2,4);
histogram(cVec,-1:0.1:1, 'FaceColor','b')
hold on
histogram(cVecBaseline,-1:0.1:1,'FaceColor',0.5*[1 1 1])
xlabel('Correlation coefficient');
ylabel('Histogram counts');
title(sprintf('S phase; %d instances',length(cVec)));
L = legend({'Performance when alone','Null model'});
L.Box = 'off'; L.Location = 'North';
L.Units = 'centimeters';
L.Position = L.Position-[0 0.1 0 0];
maxHist = 290;
text(-0.9,maxHist-3,'The predictive power of...','VerticalAlignment','Top','FontSize',9)
axis([-1 1 0 maxHist])

[cVec, cVecBaseline] = evaluatePredictivePower(d.infoV);
subplot(2,2,1);
sel = d.infoV.success(1,:)>0; %<=d.infoS.survCount(1);
rkAbd = getRank(-d.infoV.success(1,:));
rkPerf = getRank(-d.infoV.performanceExtEnv(1,:));
plot(rkPerf(~sel),rkAbd(~sel),'.','MarkerSize',15,'Color',[1 0.6 0.6])
hold on;
plot(rkPerf(sel),rkAbd(sel),'.','MarkerSize',15,'Color',[1 0 0])
xlabel('Performance rank when alone');
ylabel('Actual success rank');
title(sprintf('V-phase example; c = %.2f ',cVec(1)));
axis([0 d.infoV.totalSpeciesCount 0 d.infoV.totalSpeciesCount]);
%
text(off,dh,'(a)',opts{:})
text(off+dx,dh,'(b)',opts{:})
%
subplot(2,2,2);
histogram(cVec,-1:0.1:1,'FaceColor','r')
hold on
histogram(cVecBaseline,-1:0.1:1,'FaceColor',0.5*[1 1 1])
xlabel('Correlation coefficient');
ylabel('Histogram counts');
title(sprintf('V phase; %d instances',length(cVec)));
L = legend({'Performance when alone','Null model'});
L.Box = 'off'; L.Location = 'North';
L.Units = 'centimeters';
L.Position = L.Position-[0 0.1 0 0];
axis([-1 1 0 maxHist])
text(-0.9,maxHist-3,'The predictive power of...','VerticalAlignment','Top','FontSize',9)
% fprintf('Correlation:  %.2f +- %.2f\n', mean(cVec),std(cVec));
% fprintf('Baseline:     %.2f +- %.2f\n', mean(cVecBaseline),std(cVecBaseline));
%%
if params.saveFigs
    fname = sprintf('FigureSXX_VS_Phase_%s', params.idStr);
    saveas(gcf,[fname '.pdf'])
    saveas(gcf,[fname '.fig'])
    saveas(gcf,[fname '.png'])
end
%%

% %% Same sup figure, but using the mutual information-based metric
% histBins = 0:0.01:0.11;
% clf;
% W = 35;
% H = 10;
% set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
%     'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
%     'Units','Centimeters','Position',[1 7 W H]); 
% % Scatter plot performance rank in external environment vs abundance at equilibrium
% [~, ~, inform, informBaseline] = evaluatePredictivePower(d.infoS);
% 
% subplot(1,4,3);
% plot(d.infoS.rkAbd(1,:),d.infoS.rkExtEnv(1,:),'b.','MarkerSize',20)
% ylabel('Performance rank when alone');
% xlabel('Abundance rank at equilibrium');
% title(sprintf('S-phase example; I = %.2f ',inform(1)));
% axis([0 50 0 500]);
% 
% subplot(1,4,4);
% histogram(inform,histBins, 'FaceColor','b')
% hold on
% histogram(informBaseline,histBins,'FaceColor',0.5*[1 1 1])
% xlabel('Correlation coefficient');
% ylabel('Histogram counts');
% title(sprintf('S phase; %d instances',length(inform)));
% L = legend({'Performance when alone','Null model'});
% L.Box = 'off'; L.Location = 'North';
% L.Units = 'centimeters';
% L.Position = L.Position-[0 0.1 0 0];
% text(-0.9,34.8,'The predictive power of...','VerticalAlignment','Top','FontSize',9)
% %axis([-1 1 0 35])
% 
% [~, ~, inform, informBaseline] = evaluatePredictivePower(d.infoV);
% subplot(1,4,1);
% plot(d.infoV.rkAbd(1,:),d.infoV.rkExtEnv(1,:),'r.','MarkerSize',20)
% ylabel('Performance rank when alone');
% xlabel('Abundance rank at equilibrium');
% title(sprintf('V-phase example; I = %.2f ',inform(1)));
% axis([0 50 0 500]);
% 
% subplot(1,4,2);
% histogram(inform,histBins,'FaceColor','r')
% hold on
% histogram(informBaseline,histBins,'FaceColor',0.5*[1 1 1])
% xlabel('Correlation coefficient');
% ylabel('Histogram counts');
% title(sprintf('V phase; %d instances',length(inform)));
% L = legend({'Performance when alone','Null model'});
% L.Box = 'off'; L.Location = 'North';
% L.Units = 'centimeters';
% L.Position = L.Position-[0 0.1 0 0];
% %axis([-1 1 0 45])
% text(-0.9,44.5,'The predictive power of...','VerticalAlignment','Top','FontSize',9)
% % fprintf('Correlation:  %.2f +- %.2f\n', mean(cVec),std(cVec));
% % fprintf('Baseline:     %.2f +- %.2f\n', mean(cVecBaseline),std(cVecBaseline));
% %%
% if params.saveFigs
%     fname = sprintf('FigureSXX_VS_Phase_MutInfo_%s', params.idStr);
%     saveas(gcf,[fname '.pdf'])
%     saveas(gcf,[fname '.fig'])
% end

end

function [cVec, cVecBaseline] = evaluatePredictivePower(info)
    repCount = size(info.harvest,1);
    cVec = zeros(1,repCount);
    cVecBaseline = zeros(1,repCount);
    for t=1:repCount
        cVec(t) = corr(info.success(t,:)', info.performanceExtEnv(t,:)','type','Spearman');
        cVecBaseline(t) = corr(info.success(t,:)', -info.cpe(t,:)','type','Spearman');
    end
end

% function [cVec, cVecBaseline, inform, informBaseline] = evaluatePredictivePower(info)
%     repCount = size(info.harvest,1);
%     N = size(info.harvest,2);
%     cVec = zeros(1,repCount);
%     cVecBaseline = zeros(1,repCount);
%     inform = zeros(1,repCount);
%     informBaseline = zeros(1,repCount);
%     for t=1:repCount
%         surv = ~isnan(info.rkAbd(t,:));
%         cVec(t) = corr(info.rkAbd(t,surv)', info.rkExtEnv(t,surv)');
%         cVecBaseline(t) = corr(info.rkAbd(t,surv)', info.rkCost(t,surv)');
%         
%         % a better way to define "predictive power"
%         % We expect of order N survivors, and we have a certain ranking
%         % (either based on actual performance, or the null model based solely on cost)
%         % How predictive is one of the other? 
%         % Let's compute the mutual information between them.
%         %
%         % Note that this measure ignores the actual abundance rank, all we
%         % care about is survival / extinction.
%         inform(t) = getMutInfo(info.rkCost(t,surv), N, info.totalSpeciesCount);
%         informBaseline(t) = getMutInfo(info.rkExtEnv(t,surv), N, info.totalSpeciesCount);
%     end
% end
% 
% function info = getMutInfo(rk, N, total)
% % survived, and in the predicted top N: 
% p=zeros(2);
% p(1,1) = sum(rk <= N);
% % did not survive, although was in the predicted top N: 
% p(2,1) = N-p(1,1);
% % survived, but was not in the the predicted top N:
% p(1,2) = sum(rk > N);
% % did not survive, and was not in the the predicted top N:
% p(2,2) = total - sum(p(:));
% info = mutinfo(p);
% end
% 
% function info = mutinfo(p)
% p = p/sum(p(:));
% marg1 = sum(p,1);
% marg2 = sum(p,2);
% indep = marg2 * marg1;
% terms = p.*log2(p./indep);
% terms(isnan(terms))=0;
% info = sum(terms(:));
% end



function generate_SupFig_CorrelatedCosts(params)
% Evaluate the effect of correlated costs
% Panel A: a more complex cost model, where the survivors tend to be
% low-cost outliers and uncorrelated with the strategy
% Panel B: make an analog of Fig. 2C and observe a qualitative similarity
% Specifically, the goal is to show that a REM-like model gives a good
% approximation for the behavior of a fully non-random model, if the p and
% epsilon are adjusted correctly
%%
%p = correlatedCostParams;
    p.N = 50;
    p.sparsity = 0.5;
    p.capacity = 100*ones(p.N,1);

    % Pick a large number of organisms to generate
    alphaMax = 1000;
    p.alpha = alphaMax;
    p.P = p.alpha*p.N;
    p.seed = 1;
    p.epsilon = 0;

    % A correlated cost model:
    p.c0 = 10;
    p.Jmag = 0.02;
    p.bias = 0.015;

prop = generateOrganisms(p);
seed = 0;
prop.cost = getStructuredCosts(p,prop,seed);
prop.capacity = p.capacity;

% compute parameter values that approximate the complicated model with a
% simpler (uncorrealted) one
effectiveParamsTheory = getEffectiveParamsTheory(p);
effectiveParamsData = getEffectiveParamsData(prop);
%%
[sim, theor] = getSimAndTheory(p,prop,effectiveParamsData);
%%
W = 14;
H = 7.1;
clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
    'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
    'Units','Centimeters','Position',[14 7 W H]); 

%%% PANEL A %%%
axes('Units','Centimeters',...
    'Position',[1.5, 1.3, 5, 5]);
plotPanelA(prop,effectiveParamsTheory, effectiveParamsData);
axis([0,prop.N,0.94 1.24]);
text(-1.3,5.5,'(a)','fontsize',18,'units','centimeters');
text(5.5,5.5,'(b)','fontsize',18,'units','centimeters');
title('Cost model approximation');
set(gca,'YTick',0.9:0.1:1.2);
adjustSizes(gca,1,12);
xlabel('$|\vec\sigma_\mu|$','Interpreter','LaTeX');
ylabel('$\chi_\mu/|\vec\sigma_\mu|$','Interpreter','LaTeX');
annotation(gcf,'arrow',[0.300868055555556 0.294568452380953],...
    [0.215510172143975 0.243838028169014],'Color','r');


%%% PANEL B %%%
axes('Units','Centimeters',...
    'Position',[8.5, 1.3, 5, 5]);
%plot(sim.alpha, sim.surv');
hold on;
mn = mean(sim.surv);
sd = std(sim.surv);
upper = mn+sd;
lower = mn-sd;
patch([sim.alpha, fliplr(sim.alpha)],[upper, fliplr(lower)],'r','FaceAlpha',0.2,'EdgeColor','none');
    
plot(theor.alpha, theor.surv,'k-','LineWidth',2);
title('Survivors');
xlabel('\alpha');
box on;
adjustSizes(gca,1,12);
%%
if params.saveFigs
    fname = sprintf('FigureSXX_VS_CostApprox_%s', params.idStr);
    saveas(gcf,[fname '.pdf'])
    saveas(gcf,[fname '.fig'])
    saveas(gcf,[fname '.png'])
end

% %%
% clf;
% k = effectiveParamsData.p*prop.N;
% y = prop.cost(prop.enzCount==k)/k;
% histogram(y,'Normalization','pdf');
% a=axis;
% hold on;
% xs = 1:0.01:2;
% plot(xs,pdf('Normal',xs,effectiveParamsData.meanCost, effectiveParamsData.meanCost*effectiveParamsData.epsilon/k));
% axis(a);
% %%
% clf;
% y = prop.cost./prop.enzCount;
% histogram(y,'Normalization','pdf');
% a=axis;
% hold on;
% xs = 1:0.01:2;
% plot(xs,pdf('Normal',xs,effectiveParamsData.meanCost, effectiveParamsData.meanCost*effectiveParamsData.epsilon/k));
% axis(a);
% %%
% clf;
% histogram(prop.enzCount,'Normalization','pdf')
% set(gca,'YScale','log')
% hold on;
% xs = 1:prop.N;
% plot(xs,pdf('Binomial',xs,prop.N,effectiveParamsData.p)*effectiveParamsData.PoverS);

%%
end

function plotPanelA(prop,epTheor, epData)
%%
    mc = epData.meanCost;
    plot(prop.enzCount, prop.cost./prop.enzCount/mc,'b.')
    axis([0,prop.N,0.95,1.24])
    hold on;
    plot(epTheor.envelope.x, epTheor.envelope.y1/mc, 'b-');
    plot(epTheor.envelope.x, epTheor.envelope.y2/mc, 'b-');
    
%     envTheor = effectiveParams2envelope(epTheor, prop.P);
%     col = 0.2*[1 1 1];
%     plot(envTheor.x, envTheor.y1/mc, '--','Color',col);
%     plot(envTheor.x, envTheor.y2/mc, '--','Color',col);
    
    envData = effectiveParams2envelope(epData, prop.P);
    plot(envData.x, envData.y1/mc, 'r-','LineWidth',1);
    plot(envData.x, envData.y2/mc, 'r-','LineWidth',1);
%%
end
   
function [sim, theor] = getSimAndTheory(p,prop, ep)
    % Compute equilibria for species with a correlated cost model
    alphaList = [1:2:20, 25:5:40];
    fname = 'simData.mat';
    if exist(fname,'file')
        sim = load(fname);
        sim = sim.sim;
    else
        tryN = 50;
        %SList = p.N*(100:100:alphaMax);
        survList = zeros(tryN,length(alphaList));
        for t = 1:tryN
            fprintf('%d/%d\n',t,tryN)
            prop.cost = getStructuredCosts(p,prop,t);
            for i=1:length(alphaList)
                S = round(alphaList(i)*prop.N);
                subprop = prop;
                subprop.enzymesInSpecies = subprop.enzymesInSpecies(1:S,:);
                subprop.cost = subprop.cost(1:S);
                subprop.enzCount = subprop.enzCount(1:S);
                [~, Delta] = findEquilibrium(subprop);
                survList(t,i) = sum(abs(Delta)<1e-8);
            end
        end
        sim.alpha = alphaList;
        sim.surv = survList;
        save(fname,'sim');
    end
    %%
    % Compare with theory (using effective parameters)
    p.epsilon = ep.epsilon;
    alphaListTh = linspace(min(alphaList), max(alphaList), 100);
    p.sparsity = ep.p;
    p.R2=0;
    survTheor = zeros(1,length(alphaListTh));
    for a=1:length(alphaListTh)
        p.alpha = alphaListTh(a);
        par = solveForParams(p);
        survTheor(a) = par.survivorFraction;
    end
    theor.alpha = alphaListTh;
    theor.surv = prop.N*survTheor;
end

function p = correlatedCostParams
    p.N = 50;
    p.sparsity = 0.5;
    p.capacity = 100*ones(p.N,1);

    % Pick a large number of organisms to generate
    alphaMax = 50;
    p.alpha = alphaMax;
    p.P = p.alpha*p.N;
    p.seed = 1;
    p.epsilon = 0;

    % A correlated cost model:
    p.c0 = 5;
    p.Jmag = 0.02;
    p.bias = 0.015;
end

function cost = getStructuredCosts(costParams,prop,seed)
    rng(seed); 
    J = costParams.bias+costParams.Jmag*randn(costParams.N);
    % now set costs using matrix J
    enzCount = sum(prop.enzymesInSpecies,2);
    cost = zeros(size(enzCount));

    for s=1:size(prop.enzymesInSpecies,1)
        cost(s) = costParams.c0 + enzCount(s) + prop.enzymesInSpecies(s,:) * J * (prop.enzymesInSpecies(s,:)');
    end
end

function env = effectiveParams2envelope(ep, S)
    p = ep.p;
    N = ep.N;
    k = p*N;
    theoryX = 1:0.01:N;
    speciesCount = round(S*ep.PoverS);
    tryN = nchoosek(N,theoryX).*(p.^theoryX).*((1-p).^(N-theoryX))*speciesCount;
    theoryY1 = ep.meanCost*(1 + ep.epsilon/k*minOfGaussianSamples(tryN));
    theoryY2 = ep.meanCost*(1 - ep.epsilon/k*minOfGaussianSamples(tryN));
    sel = theoryY1<theoryY2;

    env.x = theoryX(sel);
    env.y1 = theoryY1(sel);
    env.y2 = theoryY2(sel);
end
