function [CSTR,Clusters] = RunClusterPipelineProp(Properties,Ncluster,makePlot,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS)

% INPUTS:
%
% Properties = a matrix of properties for each pipeline, where each row is
% a property and each column is a pipeline
%
% Ncluster = number of clusters to extract
%
% makePlot = set to 1 to make the plot
%
% PROCESSING_MATRIX = A matrix indicating which processing choices were
% used
%
% PROCESSING_MATRIX_LABELS = A cell wehere each cell indicates what each
% row of PROCESSING_MATRIX is
%
% OUTPUTS:
%
% CSTR = Correlation matrix generated based on the correlations between
% pipelines (note this is in the same order as the pipelines are in 
% the Properties input
%
% Clusters = The cluster number of each pipeline

n = size(Properties,2)+1;

CSTR = corr(Properties,'Type','Spearman');
[IDX,~,~,tree] = BF_ClusterReorder(CSTR);

% dendrogram_pos = [.1    0.15    0.085    0.8];
% 
% matrix_pos = [0.2019    0.15    0.6981    0.8];
% 
% cmap_pos = [0.9073    0.15    0.0171    0.8];
% 
% ypos = .05;
% legendplot = [0.2019    ypos    0.6981    0.15-ypos];

XPOS = .1319;

XPOS_OFFSET = .03;

XLENGTH = 0.7981;

YLENGTH = .825;

dendrogram_pos = [XPOS_OFFSET    0.15    0.085    YLENGTH];

matrix_pos = [XPOS    0.15    XLENGTH    YLENGTH];

cmap_pos = [XPOS+XLENGTH+.01    0.15    0.0171    YLENGTH];

ypos = .05;
legendplot = [XPOS    ypos    XLENGTH    0.15-ypos];

Clusters = cluster(tree,'maxclust',Ncluster);
cutoff = median([tree(end-(Ncluster-2),3) tree(end-(Ncluster-1),3)]);

if makePlot
figure('Position',[0 0 2407 1887]);


%ax_matrix = subplot_tight(88,8,subplot_inds1);
ax_matrix = axes('Position',matrix_pos);
%ax_matrix.Position = [0.2019    0.15    0.6981    0.85];
imagesc(CSTR(IDX,IDX))

minCorr = min(CSTR(:));
maxCorr = max(CSTR(:));

if minCorr < 0 && maxCorr > 0

caxis([-1 1])

cmap = [make_cmap('steelblue',250,30,0);flipud(make_cmap('orangered',250,30,0))];

elseif minCorr >= 0 && maxCorr >= 0
    
    caxis([0 1])
    cmap = flipud(make_cmap('orangered',250,30,0));
    
elseif minCorr <= 0 && maxCorr <= 0
      caxis([-1 0])
    cmap = make_cmap('steelblue',250,30,0);  
end

colormap(cmap)

c = colorbar;
c.Position = cmap_pos;
c.Label.String = 'Correlation between node strength';
c.FontSize = 20;
yticks(1:80)
yticklabels(IDX)
xticks([])
%yticks([])
set(gca,'FontSize',18);
ax_matrix.TickLength = [0 0];

%ax_dendrogram = subplot_tight(88,8,subplot_inds2);
ax_dendrogram = axes('Position',dendrogram_pos);

[H,~,outperm] = dendrogram(tree,0,'Orientation','left','ColorThreshold',cutoff,'Reorder',fliplr(IDX));
set(H,'LineWidth',2)
xticks([])
yticks([])

lineColours = cell2mat(get(H,'Color'));
colourList = unique(lineColours, 'rows');

%% Assign descired colours to dendrogram and clustering.
% There is definitely a smarter way to do this but this is what I came up
% with

myColours = [247,129,191; 152,78,163; 77,175,74; 55,126,184; 228,26,28;...
255,127,0;...
255,255,51;...
166,86,40;...
153,153,153]/255;

myColours = repmat(myColours,3,1);

% Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    % Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    % Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
% Apply the new colours to the chart's line objects (line by line)
for line = 1:size(H,1)
    set(H(line), 'Color', lineColours(line,:));
end
axis off

% This extracts the colours each original index was assigned

color_assigned = zeros(n-1,3);

for i=1:size(H,1)
    for j=H(i).YData(H(i).XData==0)
        color_assigned(outperm(j),:)=H(i).Color;
    end
end

% Find the cluster membership of the reordered obervations and the colour
% assigned to them

IDX_clusters = Clusters(IDX);
IDX_clusters_color = color_assigned(IDX,:);

ClusterColour = zeros(Ncluster,3);

% Draw rectangles

for i = 1:Ncluster
    IDXposition = find(IDX_clusters==i);
    
    clust_start_xy = min(IDXposition);
    clust_size = length(IDXposition);
    
    % Get the colour assigned to each cluster
    
    ClusterColour(i,:) = IDX_clusters_color(clust_start_xy,:);
    
    % Through trial and error I found these values work best. If the figure
    % size changes however these will not give the best results
    
    offset1 = 1/3;
    offset2 = 1/3;
    
    % Plot a black rectangle first to add some contrast
    
    rectangle(ax_matrix,'Position',[clust_start_xy-offset1 clust_start_xy-offset1 clust_size-offset2 clust_size-offset2],'EdgeColor','k',...
            'LineWidth',6);
        
    % Plot the coloured rectangle corresponding to the colour of that
    % cluster
    
    rectangle(ax_matrix,'Position',[clust_start_xy-offset1 clust_start_xy-offset1 clust_size-offset2 clust_size-offset2],'EdgeColor',ClusterColour(i,:),...
            'LineWidth',3.75);
    
end


%ax1 = subplot_tight(88,8,subplot_inds);

ax1 = axes('Position',legendplot);

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;

imagesc(PROCESSING_MATRIX(:,IDX))
yticks(1:length(PROCESSING_MATRIX_LABELS));
yticklabels(PROCESSING_MATRIX_LABELS);
set(gca, 'YTickLabel', PROCESSING_MATRIX_LABELS);
ytickangle(0)
xticks(1:n-1)
xtickangle(90);

xticklabels(IDX)

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

set(gca,'FontSize',18);


end

end