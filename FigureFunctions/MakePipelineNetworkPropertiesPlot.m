function MakePipelineNetworkPropertiesPlot(PARCELLATION)

% Generate a figure of some basic properties for each pipeline

% INPUT:
%
% PARCELLATION = set to 1, 2 or 3 for the 82 node, 220 node, or 380 node
% parcellation respectively

load('Pipelines_EdgeProperties_thr_0.05_inc0Edges_0.mat')

[ORDERED_INDS,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS] = FindPipelineCombinations([0 0 0 0 0 0 PARCELLATION],[7 1 6 2 4 3 5],1);

figure('Position',[0 0 2560 1440]);

subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.05], [0.1 0.01], [0.1 0.1]); 

subplot_tight(88,1,39:65);

ax = gca;

right_color = [248 184 139]./255;
left_color = [119 158 203]./255;
set(ax,'defaultAxesColorOrder',[right_color; left_color]);

FA_inds = [21:40 61:80];
NOS_inds = [1:20 41:60];

meanedgeweight_mean_NOS = nan(1,80);

meanedgeweight_mean_FA = nan(1,80);

meanedgeweight_mean_NOS(NOS_inds) = MeanEdgeWeight_mean(ORDERED_INDS(NOS_inds));

meanedgeweight_mean_FA(FA_inds) = MeanEdgeWeight_mean(ORDERED_INDS(FA_inds));

yyaxis left

bar(meanedgeweight_mean_NOS,'BarWidth', 1,'FaceColor',[119 158 203]./255,'EdgeColor','k');

xlim([0.5 80.5])

ylim([min(ylim) max(ylim)*1.1])
box on
ytickangle(90)
xticks([]);
ylabel({'Mean NOS edge weight'})
set(gca,'FontSize',16);

yyaxis right

bar(meanedgeweight_mean_FA,'BarWidth', 1,'FaceColor',[248 184 139]./255,'EdgeColor','k');

xlim([0.5 80.5])
ylim([min(ylim) max(ylim)*1.1])
box on
ytickangle(90)
xticks([]);
ylabel({'Mean FA edge weight'})
set(gca,'FontSize',16);


den_vals = density_mean(ORDERED_INDS);

if max(den_vals) > 1 
    den_vals = den_vals ./length(triu2vec(EdgeMatWeight{ORDERED_INDS(1)},1));
end


subplot_tight(88,1,12:38)
bar(den_vals,'BarWidth', 1,'FaceColor',[243 106 103]./255,'EdgeColor','k');
xticks([]);
xlim([0.5 80.5]);
ytickangle(90)
ylim([0 1.1])
box on
%ylim([0 .25])

ylabel({'Mean network density'})

set(gca,'FontSize',16);

ax1 = subplot_tight(88,1,66:77);

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;

n = size(PROCESSING_MATRIX,2)+1;

imagesc(PROCESSING_MATRIX)
yticks(1:length(PROCESSING_MATRIX_LABELS));
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

Alabel = annotation('textbox',[0, 0.8637, 0.0160, 0.0280],'String','A','EdgeColor','none','FontSize',32);

Blabel = annotation('textbox',[0, 0.5881, 0.0160, 0.0280],'String','B','EdgeColor','none','FontSize',32);

Clabel = annotation('textbox',[0, 0.3125, 0.0160, 0.0280],'String','C','EdgeColor','none','FontSize',32);
