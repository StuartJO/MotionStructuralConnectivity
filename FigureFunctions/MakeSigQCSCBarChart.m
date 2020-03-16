function MakeSigQCSCBarChart(QCSC_DATA,QCSC_PVALS,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS)

% This function makes a bar chart showing the proportion of edges in each
% pipeline that had a significant QC-SC correlation, and the proportion of
% those which were positive and negative
%
% INPUTS:
%
% QCSC_DATA = A cell array where each cell contains a vector of each edges QCSC correlation p value
%
% QCSC_PVALS = A cell array where each cell contains a vector of each edges QCSC correlation p value
%
% PROCESSING_MATRIX = A matrix indicating which processing choices were
% used
%
% PROCESSING_MATRIX_LABELS = A cell wehere each cell indicates what each
% row of PROCESSING_MATRIX is

n = size(PROCESSING_MATRIX,2)+1;

for i = 1:length(QCSC_PVALS)
    SIGEDGES = logical((QCSC_PVALS{i} < .05));
    NEGSIG = length(find(QCSC_DATA{i}(SIGEDGES) < 0));
    prop_sig(i) = length(find(SIGEDGES))/length(QCSC_PVALS{i});
    prop_sig_neg(i) = NEGSIG/length(QCSC_PVALS{i});
end

figure('Position',[0 0 2560 1440]);


subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.05], [0.1 0.01], [0.1 0.01]); 

subplot_tight(88,1,36:65)
extraParams.customSpot = '';
extraParams.add0Line = true;

for i = 1:length(QCSC_DATA)
    theColors{i} = [248 184 139]./255;
end

extraParams.theColors = theColors;

data = struct;
for i = 1:length(QCSC_DATA)
    data.(['Data',num2str(i)]) = QCSC_DATA{i};
end
    
violinplot(data,[],'ViolinAlpha',1,'EdgeColor',[248 184 139]./255,'ViolinColor',[248 184 139]./255,'BoxColor',[248 184 139]./255,'ShowData',false,'MedianColor',[248 184 139]./255,'ShowMean',true);

hold on
plot([.5 80.5],[0 0],'k--','LineWidth',1)
xlim([0.5 80.5])
ylim([-1.2 1.2])
box on
ytickangle(90)
xticks([]);
ylabel({'QC-SC correlation'})
set(gca,'FontSize',16);

subplot_tight(88,1,12:35)

PosBar = bar(0,'FaceColor',[243 106 103]./255,'EdgeColor','k');
hold on
NegBar = bar(0,'FaceColor',[119 158 203]./255,'EdgeColor','k');

lg = legend('Positive Correlation','Negative Correlation','Location','northeast','AutoUpdate','off');
title(lg,{'Proportion of significant','edges with a'}) 
lg.Box = 'off';

%delete(PosBar)
%delete(NegBar)

bar(prop_sig,'BarWidth', 1,'FaceColor',[243 106 103]./255,'EdgeColor','k');
hold on
bar(prop_sig_neg,'BarWidth', 1,'FaceColor',[119 158 203]./255,'EdgeColor','k');
xticks([]);
xlim([0.5 n-.5]);
ytickangle(90)
ylim([0 .85])
box on

ylabel({'Proportion of edges with a','significant QC-SC correlation'})

set(gca,'FontSize',16);

ax1 = subplot_tight(88,1,66:77);

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;

imagesc(PROCESSING_MATRIX)
yticks(1:length(PROCESSING_MATRIX_LABELS));
yticklabels(PROCESSING_MATRIX_LABELS);
set(gca, 'YTickLabel', PROCESSING_MATRIX_LABELS);
ytickangle(0)
xticks(1:n-1)
xtickangle(90);
hold on

for i = 0:1:n
    plot([0 n],[i-.5 i-.5],'k')
end

for i = 0:n
    plot([i-.5 i-.5],[0 n],'k')
end

xlabel('Pipeline')

colormap(ax1,[Color1; Color2; Color3])

ax = gca;

ax.TickLength = [0 0];

set(gca,'FontSize',12);

Alabel = annotation('textbox',[0.0117, 0.8637, 0.0160, 0.0280],'String','A','EdgeColor','none','FontSize',32);

Blabel = annotation('textbox',[0.0117, 0.6177, 0.0160, 0.0280],'String','B','EdgeColor','none','FontSize',32);

Clabel = annotation('textbox',[0.0117, 0.3125, 0.0160, 0.0280],'String','C','EdgeColor','none','FontSize',32);

end