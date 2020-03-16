% Note this assumes you are in the location where this script is located
PATH = pwd;

addpath(genpath(PATH))

% Set to 1 to save the figures as they are made
saveFigs = 1;

% Load in the data
load('.\analysed_data\Pipelines_EdgeProperties_thr_0.05_inc0Edges_0.mat')
load('.\analysed_data\Pipelines_QCSC_thr_0.05_RMS_abs_inc0Edges_0.mat')
load('.\analysed_data\Node_degree_strength_thr_0.05_inc0Edges_0.mat')

% Get the ordering of the data
[ORDERED_INDS_220,ORDERED_MATRIX_220,LABELS_220] = FindPipelineCombinations([0 0 0 0 0 0 2],[7 1 6 2 4 3 5],1);
[ORDERED_INDS_82,ORDERED_MATRIX_82,LABELS_82] = FindPipelineCombinations([0 0 0 0 0 0 1],[7 1 6 2 4 3 5],1);
[ORDERED_INDS_380,ORDERED_MATRIX_380,LABELS_380] = FindPipelineCombinations([0 0 0 0 0 0 3],[7 1 6 2 4 3 5],1);

% Figure 2
MakePipelineNetworkPropertiesPlot(2)

if saveFigs == 1
    print('Fig1.tif','-dtiff','-r300')
end

% Figure S2
MakePipelineNetworkPropertiesPlot(1)

if saveFigs == 1
    print('FigS2.tif','-dtiff','-r300')
end

% Figure S3
MakePipelineNetworkPropertiesPlot(3)

if saveFigs == 1
    print('FigS3.tif','-dtiff','-r300')
end

% Figure 3

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('Fig3.tif','-dtiff','-r300')
end

% Figure S4

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_82),QCSC_PVALS(ORDERED_INDS_82),ORDERED_MATRIX_82,LABELS_82)

if saveFigs == 1
    print('FigS4.tif','-dtiff','-r300')
end

% Figure S5

MakeSigQCSCBarChart(QCSC(ORDERED_INDS_380),QCSC_PVALS(ORDERED_INDS_380),ORDERED_MATRIX_380,LABELS_380)

if saveFigs == 1
    print('FigS5.tif','-dtiff','-r300')
end


% Figure 4
QCSCvsEdgeConLenVar(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),[11 14 34 46 80])

if saveFigs == 1
    print('Fig4.tif','-dtiff','-r300')
end

% Figure S15

QCSCvsEdgeConLenWei(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeight(ORDERED_INDS_220),[11 14 34 46 80])

if saveFigs == 1
    print('FigS15.tif','-dtiff','-r300')
end

% Figure S7

QCSCperThreshold(QCSC_PVALS(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('FigS7.tif','-dtiff','-r300')
end

% Figure 5

QCStrengthperThreshold(PropNodeStr_con_sig(ORDERED_INDS_220,:),PropNodeStr_var_sig(ORDERED_INDS_220,:),ORDERED_MATRIX_220,LABELS_220)

if saveFigs == 1
    print('Fig5.tif','-dtiff','-r300')
end

% Figure 6 subplots. These were subsequently combined in PowerPoint to make
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

% Figure 7: clustering

mean_strength_rank = zeros(220,80);

for i = 1:length(ORDERED_INDS_220)
    
    STR = STRcon{ORDERED_INDS_220(i),7};
    mean_strength_rank(:,i) = mean(tiedrank(STR'),2)';

end

[CTSR,clust] = RunClusterPipelineProp(mean_strength_rank,4,1,ORDERED_MATRIX_220,LABELS_220);

if saveFigs == 1
    print('Fig7.tif','-dtiff','-r300')
end

[CorrVec,CorrIND] = triu2vec(CTSR,1);

[minCorr,minCorrIND] = min(CorrVec);

[Pipe1,Pipe2] = ind2sub(size(CTSR),CorrIND(minCorrIND));

% Figure 8, panel A

data = mean_strength_rank(:,Pipe1);

PlotBrainNodalProperty(data,['Pipeline ',num2str(Pipe1)],'A',[0 220],3,0,[0 220],0)

if saveFigs == 1
    print('Fig8_A.tif','-dtiff','-r300')
end

% Figure 8, panel B

data = mean_strength_rank(:,Pipe2);

PlotBrainNodalProperty(data,['Pipeline ',num2str(Pipe2)],'B',[0 220],3,0,[0 220],0)

if saveFigs == 1
    print('Fig8_B.tif','-dtiff','-r300')
end

strdatasigcvals = zeros(220,80);

for i = 1:80
    C = QCSTR_con{ORDERED_INDS_220(i),7};
    P = QCSTR_pval_con{ORDERED_INDS_220(i),7};
    strdatasigcvals(:,i) = C.*(P<.05);
end

NodePropSig = mean(strdatasigcvals > 0,2);

% Figure S17

PlotBrainNodalProperty(NodePropSig,'','',[0 1],2,'Proportion significant',[0 1])

if saveFigs == 1
    print('FigS17.tif','-dtiff','-r300')
end

% Figure S1

PlotMotionProperties

if saveFigs == 1
    print('FigS1.tif','-dtiff','-r300')
end


% Figure S6

CompareEDDYversions

if saveFigs == 1
    print('FigS6.tif','-dtiff','-r300')
end


% Figures S8 to S14

for i = 2:length(MOTIONNAMES)
    load(['.\analysed_data\Pipelines_QCSC_thr_0.05_',MOTIONNAMES{i},'_inc0Edges_0.mat'])
    MakeSigQCSCBarChart(QCSC(ORDERED_INDS_220),QCSC_PVALS(ORDERED_INDS_220),ORDERED_MATRIX_220,LABELS_220)
    if saveFigs == 1
        print(['FigS',num2str(i+6),'.tif'],'-dtiff','-r300')
    end

end


% Figure S16

load('.\analysed_data\Pipelines_EdgeProperties_thr_0_inc0Edges_1.mat')
load('.\analysed_data\Pipelines_QCSC_thr_0_RMS_abs_inc0Edges_1.mat')

QCSCvsEdgeConLenVar(QCSC(ORDERED_INDS_220),EdgeLength(ORDERED_INDS_220),EdgeConsistency(ORDERED_INDS_220),EdgeWeightVariability(ORDERED_INDS_220),[11 14 34 46 80])

if saveFigs == 1
    print('FigS16.tif','-dtiff','-r300')
end

