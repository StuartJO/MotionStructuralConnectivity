function [prop_sig,ax1,ax2,ax3] = MakeSigQCSCBarChart(QCSC_DATA,QCSC_PVALS,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS)

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

%figure('Position',[0 0 2560 1440]);
%figure('Position',[0 0 2560 1052]);

figure('Position',[0 0 2560 1440]);

PanelAPos = [0.146  0.6185  0.8476  0.3650];
PanelBPos = [0.146  0.2495  0.8476  0.3650];
PanelCPos = [0.146  0.0808  0.8476  0.1647];

% PanelAPos = [0.15  0.6185  0.8476  0.3650];
% PanelBPos = [0.15  0.2495  0.8476  0.3650];
% PanelCPos = [0.15  0.0808  0.8476  0.1647];

LblYOffset = .0075;

FONTSIZE = 22;

ALblPos = [0, (PanelAPos(2)+PanelAPos(4)-LblYOffset), 0.0160, 0.0280];
BLblPos = [0, (PanelBPos(2)+PanelBPos(4)-LblYOffset), 0.0160, 0.0280];
CLblPos = [0, (PanelCPos(2)+PanelCPos(4)-LblYOffset), 0.0160, 0.0280];

%ax2 = axes('Position',[0.1424  0.3339  0.8476  0.2696]);
ax2 = axes('Position',PanelBPos);

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
set(gca,'FontSize',FONTSIZE);

%ax1 = subplot_tight(88,1,12:38);
%ax1 = axes('Position',[0.1424  0.6085  0.8476  0.2696]);
ax1 = axes('Position',PanelAPos);

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
%ylim([0 .85])
ylim([0 max(prop_sig)+.1])
box on

ylabel({'Proportion of edges with a','significant QC-SC correlation'})

set(gca,'FontSize',FONTSIZE);

%ax3 = subplot_tight(88,1,66:77);
%ax3 = axes('Position',[0.1424  0.2119  0.8476  0.1170]);
ax3 = axes('Position',PanelCPos);

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

colormap(ax3,[Color1; Color2; Color3])

ax = gca;

ax.TickLength = [0 0];

set(gca,'FontSize',FONTSIZE);

Alabel = annotation('textbox',ALblPos,'String','A','EdgeColor','none','FontSize',32);

Blabel = annotation('textbox',BLblPos,'String','B','EdgeColor','none','FontSize',32);

Clabel = annotation('textbox',CLblPos,'String','C','EdgeColor','none','FontSize',32);

end