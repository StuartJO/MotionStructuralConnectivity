%% Figure generation script

% Note that all the figures are generated at a specifc size. I generated
% these on a windows PC at 4K screen resolution with scaling set to 100% in
% windows settings. If the scaling is set differently that will effect how
% the images look (MATLAB seems to determine screen resolution based on the
% scaling of the default resolution). So if aiming to reproduce these plots
% or adapt yourself, keep that in mind


% Note this assumes you are in the location where this script is located
PATH = pwd;

%addpath(genpath(PATH))

% Set to 1 to save the figures as they are made
saveFigs = 1;

% Load in the data
load('./analysed_data/Pipelines_EdgeProperties_thr_0.05_inc0Edges_0.mat')
load('./analysed_data/Pipelines_QCSC_thr_0.05_ABSall_inc0Edges_0.mat')
load('./analysed_data/Node_degree_strength_thr_0.05_inc0Edges_0.mat')

load('COMBINATIONS_MATRIX.mat')

% Get the ordering of the data
[ORDERED_INDS_220,ORDERED_MATRIX_220,LABELS_220] = FindPipelineCombinations([0 0 0 0 0 0 2],[7 1 6 2 4 3 5],1);
[ORDERED_INDS_82,ORDERED_MATRIX_82,LABELS_82] = FindPipelineCombinations([0 0 0 0 0 0 1],[7 1 6 2 4 3 5],1);
[ORDERED_INDS_380,ORDERED_MATRIX_380,LABELS_380] = FindPipelineCombinations([0 0 0 0 0 0 3],[7 1 6 2 4 3 5],1);

% Figure 2
MakePipelineNetworkPropertiesPlot(2)

if saveFigs == 1
    print('FigS2.tif','-dtiff','-r300')
    %print('smallFig2.png','-dpng')
end

% Figure S2
MakePipelineNetworkPropertiesPlot(1)

if saveFigs == 1
    print('FigS1.tif','-dtiff','-r300')
end

% Figure S3
MakePipelineNetworkPropertiesPlot(3)

if saveFigs == 1
    print('FigS3.tif','-dtiff','-r300')
end

% Figure 2

N220_nsig = MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('Fig2.tif','-dtiff','-r300')
    print('smallFig2.tif','-dtiff')
end

N220_nsig_fdr = MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS_FDR(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('FigS5.tif','-dtiff','-r300')
    print('smallFigS5.tif','-dtiff')
end

% Figure S4

N82_nsig = MakeSigQCSCBarChart(QCSC(ORDERED_INDS_82),QCSC_PVALS(ORDERED_INDS_82),ORDERED_MATRIX_82,LABELS_82);

if saveFigs == 1
    print('FigS6.tif','-dtiff','-r300')
end

% Figure S5

N380_nsig = MakeSigQCSCBarChart(QCSC(ORDERED_INDS_380),QCSC_PVALS(ORDERED_INDS_380),ORDERED_MATRIX_380,LABELS_380);

if saveFigs == 1
    print('FigS7.tif','-dtiff','-r300')
end


% Figure 3
Inds2Use = [11 34 46 66];

Pipelines2Use = ORDERED_INDS_220(Inds2Use);

for i = 1:4
    IDX = sub2ind(size(PROCESSING_STEPS_MATRIX),1:7,COMBINATIONS(Pipelines2Use(i),:));

S = strjoin(PROCESSING_STEPS_MATRIX(IDX([6 2])),'+');
PipelineLabels{i} = [S,' (pipeline ',num2str(Inds2Use(i)),')'];

end

QCSCvsEdgeConLenVar(QCSC(Pipelines2Use),EdgeLength(Pipelines2Use),EdgeConsistency(Pipelines2Use),EdgeWeightVariability(Pipelines2Use),PipelineLabels)

if saveFigs == 1
    print('Fig3.tif','-dtiff','-r300')
    print('smallFig3.tif','-dtiff')
end

% Figure S15

QCSCvsEdgeConLenWei(QCSC(Pipelines2Use),EdgeLength(Pipelines2Use),EdgeConsistency(Pipelines2Use),EdgeWeight(Pipelines2Use),PipelineLabels)

if saveFigs == 1
    print('FigS19.tif','-dtiff','-r300')
end

% Figure S7

QCSCperThreshold(QCSC_PVALS(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('FigS10.tif','-dtiff','-r300')
end

% Figure 4

QCStrengthperThreshold(PropNodeStr_con_sig(ORDERED_INDS_220,:),PropNodeStr_var_sig(ORDERED_INDS_220,:),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('Fig4.tif','-dtiff','-r300')
    print('smallFig4.tif','-dtiff')
end

% Figure 6 subplots. These were subsequently combined in PowerPoint to make
% the image seen in the paper (I couldn't get the figures to be
% appriopriately sized through MATLAB)

inds2use = [11 14 40 51 54 80];

Pipelines2Use = ORDERED_INDS_220(inds2use);

PlotLabelLetters = {'A','B','C','D','E','F'};
climits = [-1 1];
PlotColorBar = 0;
Ignore0vals = 1;

for i = 1:length(inds2use)

    C = QCSTR_con{Pipelines2Use(i),7};
    P = QCSTR_pval_con{Pipelines2Use(i),7};
    
    
    IDX = sub2ind(size(PROCESSING_STEPS_MATRIX),1:7,COMBINATIONS(Pipelines2Use(i),:));

    %plotName = ['Pipeline ',num2str(inds2use(i))];
    
    S = strjoin(PROCESSING_STEPS_MATRIX(IDX([1 6 2 4 3 5])),', ');
    
    
    %plotName = {['Pipeline ',num2str(inds2use(i)),':',S{1},',',S{6},',',S{2}],[',',S{4},',',S{3},',',S{5}]};
plotName = ['Pipeline ',num2str(inds2use(i)),': ',S];

strsigcvals = C.*(P<.05);

T = PlotBrainNodalProperty(strsigcvals,plotName,PlotLabelLetters{i},climits,1,PlotColorBar,[],Ignore0vals);

if saveFigs == 1
    print(['Fig5_',PlotLabelLetters{i},'.tif'],'-dtiff','-r300')
end

end

% Figure 6: clustering

mean_strength_rank = zeros(220,80);

for i = 1:length(ORDERED_INDS_220)
    
    STR = STRcon{ORDERED_INDS_220(i),7};
    mean_strength_rank(:,i) = mean(tiedrank(STR'),2)';

end

[CSTR,clust] = RunClusterPipelineProp(mean_strength_rank,6,1,ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('Fig6.tif','-dtiff','-r300')
    print('smallFig6.png','-dpng')
end

for i = 1:max(clust)
    INDS = clust==i;
    clust_qcstr = triu2vec(CSTR(INDS,INDS),1);
    clust_min(i) = min(clust_qcstr);
    clust_max(i) = max(clust_qcstr);
    Csize(i) = sum(INDS);
    clust_pipes{i} = find(clust==i);
end

% Find the pipelines with the worst strength correlation

[CorrVec,CorrIND] = triu2vec(CSTR,1);

[minCorr,minCorrIND] = min(CorrVec);

[Pipe1,Pipe2] = ind2sub(size(CSTR),CorrIND(minCorrIND));

% Figure 8, panel A

data = mean_strength_rank(:,Pipe1);

    IDX = sub2ind(size(PROCESSING_STEPS_MATRIX),1:7,COMBINATIONS(ORDERED_INDS_220(Pipe1),:));

    %plotName = ['Pipeline ',num2str(inds2use(i))];
    
    S = strjoin(PROCESSING_STEPS_MATRIX(IDX([1 6 2 4 3 5])),', ');
    
    
    %plotName = {['Pipeline ',num2str(inds2use(i)),':',S{1},',',S{6},',',S{2}],[',',S{4},',',S{3},',',S{5}]};
plotName = ['Pipeline ',num2str(Pipe1),': ',S];

PlotBrainNodalProperty(data,plotName,'A',[0 220],3,0,[0 220],0)

if saveFigs == 1
    print('Fig7_A.tif','-dtiff','-r300')
end

% Figure 8, panel B

data = mean_strength_rank(:,Pipe2);

IDX = sub2ind(size(PROCESSING_STEPS_MATRIX),1:7,COMBINATIONS(ORDERED_INDS_220(Pipe2),:));

    %plotName = ['Pipeline ',num2str(inds2use(i))];
    
    S = strjoin(PROCESSING_STEPS_MATRIX(IDX([1 6 2 4 3 5])),', ');
    plotName = ['Pipeline ',num2str(Pipe2),': ',S];

PlotBrainNodalProperty(data,plotName,'B',[0 220],3,0,[0 220],0)

if saveFigs == 1
    print('Fig7_B.tif','-dtiff','-r300')
end

strdatasigcvals = zeros(220,80);

for i = 1:80
    C = QCSTR_con{ORDERED_INDS_220(i),7};
    P = QCSTR_pval_con{ORDERED_INDS_220(i),7};
    strdatasigcvals(:,i) = C.*(P<.05);
end

NodePropSig = mean(strdatasigcvals > 0,2);

% Figure S21

PlotBrainNodalProperty(NodePropSig,'Proportion significant','',[0 1],2,1,[0 1])

if saveFigs == 1
    print('FigS21.tif','-dtiff','-r300')
end

% Figure S1

PlotMotionProperties

if saveFigs == 1
    print('FigS4.tif','-dtiff','-r300')
end


% Figure S8

CompareEDDYversions

if saveFigs == 1
    print('FigS8.tif','-dtiff','-r300')
end


% Figures S11 to S18

load('MOTION_DATA.mat', 'MOTIONNAMES')
for i = 2:9
    load(['Pipelines_QCSC_thr_0.05_',MOTIONNAMES{i},'_inc0Edges_0.mat'])
    MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220);
    if saveFigs == 1
        print(['FigS',num2str(i+9),'.tif'],'-dtiff','-r300')
    end

end


% Figure S20

load('./analysed_data/Pipelines_EdgeProperties_thr_0_inc0Edges_1.mat')
load('./analysed_data/Pipelines_QCSC_thr_0_ABSall_inc0Edges_1.mat')

Inds2Use = [11 34 46 66];

Pipelines2Use = ORDERED_INDS_220(Inds2Use);

for i = 1:4
    IDX = sub2ind(size(PROCESSING_STEPS_MATRIX),1:7,COMBINATIONS(Pipelines2Use(i),:));

S = strjoin(PROCESSING_STEPS_MATRIX(IDX([6 2])),'+');
PipelineLabels{i} = [S,' (pipeline ',num2str(Inds2Use(i)),')'];

end

QCSCvsEdgeConLenVar(QCSC(Pipelines2Use),EdgeLength(Pipelines2Use),EdgeConsistency(Pipelines2Use),EdgeWeightVariability(Pipelines2Use),PipelineLabels);

if saveFigs == 1
    print('FigS20.tif','-dtiff','-r300')
end

% Figure S9

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('FigS9.tif','-dtiff','-r300')
end

   
 % For some figures I combined the images seperately in another program
 % (PowerPoint of all things!!) and the following just makes some large
 % colourbars I could use because otherwise I would have either colourbars
 % for every plot which is just silly or b) a comically small colourbar for
 % one plot.
    
figure('Position',[1 41 500 1323])
c = colorbar('Position',[0.5    0.1109    0.1    0.8135]);
axis off
c.Label.String = 'QC-Strength correlation';
c.FontSize = 28;
c.Label.FontSize = 34;
caxis([-1 1])
set(c, 'xlim', [-1 1])
cmap = [make_cmap('steelblue',250,30,0);flipud(make_cmap('orangered',250,30,0))];
colormap(cmap)
c.LineWidth = 2;
print('CorrColorBarVert.tif','-dtiff','-r300')

figure('Position',[1 41 500 1323])
c = colorbar('Position',[0.5    0.1109    0.1    0.8135]);
axis off
c.Label.String = 'Node strength rank';
c.FontSize = 28;
c.Label.FontSize = 34;
caxis([0 220])
set(c, 'xlim', [0 220])
cmap = viridisplus(500);
colormap(cmap)
c.LineWidth = 2;
print('StrColorBarVert.tif','-dtiff','-r300')


% Calculate some statistics used in the paper

load('MOTION_DATA.mat','motion_data','MOTIONNAMES')

% Get the means and standard deviations of the various motion measures
for i = 3:11 
    motionM_SD(i-2,1) = mean(motion_data{i}); 
    motionM_SD(i-2,2) = std(motion_data{i}); 
end

% Get the mean density of FACT and iFOD2 pipelines

[FACT_PIPELINES,~,~] = FindPipelineCombinations([0 1 0 0 0 0 0],[7 1 6 2 4 3 5],1);

[iFOD2_PIPELINES,~,~] = FindPipelineCombinations([0 2 0 0 0 0 0],[7 1 6 2 4 3 5],1);

FACT_den_mean = mean(density_mean(FACT_PIPELINES));

iFOD2_den_mean = mean(density_mean(iFOD2_PIPELINES));

% Get the mean, max, and min QC-strength correlations for a 50% consistency threshold

for i = 1:length(ORDERED_INDS_220)
    Cmean_STR220N(i) = nanmean(QCSTR_con{ORDERED_INDS_220(i),7});
    
    Cmax_STR220N(i) = max(QCSTR_con{ORDERED_INDS_220(i),7});
    Cmin_STR220N(i) = min(QCSTR_con{ORDERED_INDS_220(i),7});
end
