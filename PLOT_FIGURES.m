% Note this assumes you are in the location where this script is located
PATH = pwd;

addpath(genpath(PATH))

% Set to 1 to save the figures as they are made
saveFigs = 1;

%% Load in the data
load('Pipelines_EdgeProperties_thr_0.05_inc0Edges_0.mat')
load('Pipelines_QCSC_thr_0.05_ABSall_inc0Edges_0.mat')
load('Node_degree_strength_thr_0.05_inc0Edges_0.mat')

%% Get the ordering of the data
% Pull out only the pipelines that used a specific parcellation
INDS220 = [0 0 0 0 0 0 2];
INDS82 = [0 0 0 0 0 0 1];
INDS380 = [0 0 0 0 0 0 3];

% Pull out the pipelines in an assigned order, namely by parcellation, EDDY type, Edge weight, 
% tractography algorithm, spatial constraint, seeding algorithm and tract reweighting (in that order)

PIPE_ORDERING = [7 1 6 2 4 3 5];

[ORDERED_INDS_220,ORDERED_MATRIX_220,LABELS_220] = FindPipelineCombinations(INDS220,PIPE_ORDERING,1);
[ORDERED_INDS_82,ORDERED_MATRIX_82,LABELS_82] = FindPipelineCombinations(INDS82,PIPE_ORDERING,1);
[ORDERED_INDS_380,ORDERED_MATRIX_380,LABELS_380] = FindPipelineCombinations(INDS380,PIPE_ORDERING,1);

%% Figure 2

% Note MakePipelineNetworkPropertiesPlot(2) calculates properties for pipeline using the 220 node
% parcellation

MakePipelineNetworkPropertiesPlot(2)

if saveFigs == 1
    print('Fig2.tif','-dtiff','-r300')
    print('smallFig2.png','-dpng')
end

%% Figure S2
MakePipelineNetworkPropertiesPlot(1)

if saveFigs == 1
    print('FigS2.tif','-dtiff','-r300')
end

%% Figure S3
MakePipelineNetworkPropertiesPlot(3)

if saveFigs == 1
    print('FigS3.tif','-dtiff','-r300')
end

%% Figure 3
% From the QCSC variable, only the values associated with the 220 nodes pipeline are extracted

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('Fig3.tif','-dtiff','-r300')
    print('smallFig3.tif','-dtiff')
end

%% Figure S4

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_82),QCSC_PVALS(ORDERED_INDS_82),ORDERED_MATRIX_82,LABELS_82)

if saveFigs == 1
    print('FigS4.tif','-dtiff','-r300')
end

%% Figure S5

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_380),QCSC_PVALS(ORDERED_INDS_380),ORDERED_MATRIX_380,LABELS_380)

if saveFigs == 1
    print('FigS5.tif','-dtiff','-r300')
end


%% Figure 4 & Figure S15

% These pipelines were selected as they are the most prototypical
PROTO_PIPES = [11 14 34 46 80];

QCSCvsEdgeConLenVar(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),PROTO_PIPES)

if saveFigs == 1
    print('Fig4.tif','-dtiff','-r300')
    print('smallFig4.tif','-dtiff')
end

QCSCvsEdgeConLenWei(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeight(ORDERED_INDS_220),PROTO_PIPES)

if saveFigs == 1
    print('FigS15.tif','-dtiff','-r300')
end

%% Figure S7

QCSCperThreshold(QCSC_PVALS(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('FigS7.tif','-dtiff','-r300')
end

%% Figure 5

QCStrengthperThreshold(PropNodeStr_con_sig(ORDERED_INDS_220,:),PropNodeStr_var_sig(ORDERED_INDS_220,:),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('Fig5.tif','-dtiff','-r300')
    print('smallFig5.tif','-dtiff')
end

%% Figure 6 subplots. 
% These were subsequently combined in PowerPoint to make
% the image seen in the paper (I couldn't get the figures to be
% appriopriately sized through MATLAB)

inds2use = [11 14 40 51 54 80];
PlotLabelLetters = {'A','B','C','D','E','F'};
climits = [-1 1];
PlotColorBar = 0;
Ignore0vals = 1;

for i = 1:6

    C = QCSTR_con{ORDERED_INDS_220(inds2use(i)),7};
    P = QCSTR_pval_con{ORDERED_INDS_220(inds2use(i)),7};
    
    plotName = ['Pipeline ',num2str(inds2use(i))];

    strsigcvals = C.*(P<.05);

    PlotBrainNodalProperty(strsigcvals,plotName,PlotLabelLetters{i},climits,1,PlotColorBar,[],Ignore0vals)

    if saveFigs == 1
        print(['Fig6_',PlotLabelLetters{i},'.tif'],'-dtiff','-r300')
    end

end

%% Figure 7

mean_strength_rank = zeros(220,80);

for i = 1:length(ORDERED_INDS_220)
    
    STR = STRcon{ORDERED_INDS_220(i),7};
    mean_strength_rank(:,i) = mean(tiedrank(STR'),2)';

end

% CTSR is the correlation matrix (it is not ordered by optimal clusters)
% clust is the cluster assignment of the input pipelines

% Visually from the matrix for this example it looks like four is a good fit.
clusters2extract = 4;
makePlot = 1;

[CTSR,clust] = RunClusterPipelineProp(mean_strength_rank,clusters2extract,makePlot,ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('Fig7.tif','-dtiff','-r300')
    print('smallFig7.png','-dpng')
end

[CorrVec,CorrIND] = triu2vec(CTSR,1);

[minCorr,minCorrIND] = min(CorrVec);

[Pipe1,Pipe2] = ind2sub(size(CTSR),CorrIND(minCorrIND));

%% Figure 8

color_range = [0 220];
ColorMap2Use = 3;
PlotColorBar = 0;
Color0s = 0;

% Panel A

data = mean_strength_rank(:,Pipe1);

PlotBrainNodalProperty(data,['Pipeline ',num2str(Pipe1)],'A',color_range,ColorMap2Use,PlotColorBar,color_range,Color0s)

if saveFigs == 1
    print('Fig8_A.tif','-dtiff','-r300')
end

% Panel B

data = mean_strength_rank(:,Pipe2);

PlotBrainNodalProperty(data,['Pipeline ',num2str(Pipe2)],'B',color_range,ColorMap2Use,PlotColorBar,color_range,Color0s)

if saveFigs == 1
    print('Fig8_B.tif','-dtiff','-r300')
end

%% Figure S17

% First pull out how many nodes have a significant QCStrength correlation

strdatasigcvals = zeros(220,80);

for i = 1:80
    C = QCSTR_con{ORDERED_INDS_220(i),7};
    P = QCSTR_pval_con{ORDERED_INDS_220(i),7};
    strdatasigcvals(:,i) = C.*(P<.05);
end

NodePropSig = mean(strdatasigcvals > 0,2);

% Now plot those values

color_range = [0 1];
ColorMap2Use = 2;
PlotColorBar = 1;
Color0s = 0;

PlotBrainNodalProperty(NodePropSig,'Proportion significant','',color_range,ColorMap2Use,PlotColorBar,color_range,Color0s)

if saveFigs == 1
    print('FigS17.tif','-dtiff','-r300')
end

%% Figure S1

PlotMotionProperties

if saveFigs == 1
    print('FigS1.tif','-dtiff','-r300')
end


%% Figure S6

CompareEDDYversions

if saveFigs == 1
    print('FigS6.tif','-dtiff','-r300')
end


%% Figures S9 to S14

load('MOTION_DATA.mat', 'MOTIONNAMES')
for i = 2:length(MOTIONNAMES)
    load(['Pipelines_QCSC_thr_0.05_',MOTIONNAMES{i},'_inc0Edges_0.mat'])
    MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)
    if saveFigs == 1
        print(['FigS',num2str(i+7),'.tif'],'-dtiff','-r300')
    end

end

%% Figure S16

load('./analysed_data/Pipelines_EdgeProperties_thr_0_inc0Edges_1.mat')
load('./analysed_data/Pipelines_QCSC_thr_0_ABSall_inc0Edges_1.mat')

QCSCvsEdgeConLenVar(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),[11 14 34 46 80])

if saveFigs == 1
    print('FigS16.tif','-dtiff','-r300')
end

%% Figure S8

    MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)
    if saveFigs == 1
        print('FigS8.tif','-dtiff','-r300')
    end
