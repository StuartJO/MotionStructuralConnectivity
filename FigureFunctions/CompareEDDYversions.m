function CompareEDDYversions

% This function compares the different version of EDDY (which we
% imaginatively call, EDDY1, EDDY1.5, EDDY2).

load('Pipelines_EDDY1.5_QCSC_thr_0.05_RMS_abs_inc0Edges_0.mat', 'COMBINATIONS','QCSC','QCSC_PVALS')

COMBINATIONS_EDDY = COMBINATIONS;
QCSC_eddy1half = QCSC;
QCSC_PVALS_eddy1half = QCSC_PVALS;

clear QCSC QCSC_PVALS

load('Pipelines_QCSC_thr_0.05_RMS_abs_inc0Edges_0.mat', 'COMBINATIONS','QCSC','QCSC_PVALS')

[~,IDX] = ismember(COMBINATIONS(:,2:7),COMBINATIONS_EDDY(:,2:7),'rows');

INDS = find(IDX>0);

QCSC = [QCSC(INDS); QCSC_eddy1half];
QCSC_PVALS = [QCSC_PVALS(INDS); QCSC_PVALS_eddy1half];

EDDY_COMBINATIONS = [COMBINATIONS(INDS,:); COMBINATIONS_EDDY];

% EDDY1.5 in the COMBINATIONS matrix is numbered as 3, which would result
% in things being sorted around the wrong way in my normal approach, so I
% temporarily make it 1.5, and then after the sorting I change it make so
% colour maps are correctly applied

EDDY_COMBINATIONS2 = changem(EDDY_COMBINATIONS,1.5,3);

[ORDERED_INDS,EDDYTYPES_PIPELINE_MATRIX,LABELS] = FindPipelineCombinationsEDDYCompare(EDDY_COMBINATIONS2,[0 0 0 0 0 0 2],[7 1 6 2 4 3 5],1);

EDDYTYPES_PIPELINE_MATRIX(EDDYTYPES_PIPELINE_MATRIX==1.5) = 3;

QCSC_ordered = QCSC(ORDERED_INDS);
QCSC_PVALS_ordered = QCSC_PVALS(ORDERED_INDS);

for i = 1:length(QCSC_PVALS_ordered)
    SIGEDGES = logical((QCSC_PVALS_ordered{i} < .05));
    NEGSIG = length(find(QCSC_ordered{i}(SIGEDGES) < 0));
    TOTAL_SIG(i) = length(find(SIGEDGES))/length(QCSC_PVALS_ordered{i});
    NEG_SIG(i) = NEGSIG/length(QCSC_PVALS_ordered{i});
end

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;

figure('Position',[0 0 1488 555])

subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.05], [0.1 0.01], [0.1 0.01]); 
ax1 = subplot_tight(3,8,18:24);

imagesc(EDDYTYPES_PIPELINE_MATRIX)
yticks(1:length(LABELS));
set(gca, 'YTickLabel', LABELS);
ytickangle(0)

hold on
for i = 0:1:length(EDDYTYPES_PIPELINE_MATRIX)+1
    plot([0 81],[i-.5 i-.5],'k')
end

for i = 0:length(EDDYTYPES_PIPELINE_MATRIX)+1
    plot([i-.5 i-.5],[0 length(EDDYTYPES_PIPELINE_MATRIX)+1],'k')
end
%axis image

xlabel('Pipeline')

colormap(ax1,[Color1; Color2; Color3])

ax = gca;

ax.TickLength = [0 0];

set(gca,'FontSize',12);

xticks([])

ax2 = subplot_tight(3,8,[2:8 10:16]);

bar(TOTAL_SIG,'BarWidth', 1,'FaceColor',[243 106 103]./255,'EdgeColor','k');
hold on
bar(NEG_SIG,'BarWidth', 1,'FaceColor',[119 158 203]./255,'EdgeColor','k');
xticks([]);
xlim([0.5 length(EDDYTYPES_PIPELINE_MATRIX)+.5]);
ytickangle(90)
ylim([0 .45])
box on

ylabel({'Proportion of edges with a','significant QC-SC correlation'})

set(gca,'FontSize',16);

