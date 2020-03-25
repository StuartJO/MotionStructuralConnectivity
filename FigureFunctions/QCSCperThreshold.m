function QCSCperThreshold(QCSC_PVALS,EdgeConsistency,EdgeVariability,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS)

% INPUTS:
%
% QCSC_PVALS = A cell array where each cell contains a vector of each edges QCSC correlation p value
%
% EdgeConsistency = A cell array where each cell contains a vector of each edges consistency
%
% EdgeVariability = A cell array where each cell contains a vector of each
% edges weight variability
%
% PROCESSING_MATRIX = A matrix indicating which processing choices were
% used
%
% PROCESSING_MATRIX_LABELS = A cell wehere each cell indicates what each
% row of PROCESSING_MATRIX is

n = length(QCSC_PVALS)+1;

thresh = [.05 .1:.1:.9];
threshvar = .1:.1:1;

for j = 1:length(thresh)

for i = 1:length(QCSC_PVALS)
   
   edges2use = EdgeConsistency{i} > thresh(j);
   
   ptile = prctile(EdgeVariability{i},threshvar(j)*100);
   
   edges2use2 = EdgeVariability{i} < ptile;
   
   prop_sig_thr(i,j) = length(find(QCSC_PVALS{i}(edges2use) < .05))/length(QCSC_PVALS{i}(edges2use));

   prop_sig_thr2(i,j) = length(find(QCSC_PVALS{i}(edges2use2) < .05))/length(QCSC_PVALS{i}(edges2use2));

end 
    
end

[X,Y] = meshgrid(1:size(prop_sig_thr,1),1:size(prop_sig_thr,2));

X = X';
Y = Y';

%scatter(X(:),prop_sig_thr(:),50,Y(:),'filled')

figure('units','pixels','outerposition',[0 0 2560 1440])


subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.05], [0.1 0.01], [0.1 0.01]); 
subplot_tight(88,1,12:38)

for i = 1:length(thresh)
idx = find(Y==i);
scatter(X(idx),prop_sig_thr(idx),50,Y(idx),'filled');
hold on
end
colormap(flipud(make_cmap('red',10)))
ylim([0 1.1])
box on
lgd1 = legend({'= 5%','= 10%','= 20%','= 30%','= 40%','= 50%','= 60%','= 70%','= 80%','= 90%'},'Orientation','horizontal');
title(lgd1,'Edge consistency-based threshold: Percentage')
xlim([0.5 80.5])
ytickangle(90)
xticks([]);
ylabel({'Proportion of edges with a','significant QC-SC correlation'})
set(gca,'FontSize',16);

subplot_tight(88,1,39:65)
for i = length(threshvar):-1:1
idx = find(Y==i);
scatter(X(idx),prop_sig_thr2(idx),50,11-Y(idx),'filled');
hold on
end
colormap(flipud(make_cmap('red',10)))
ylim([0 1.1])
box on
lgd2 = legend(fliplr({'= 10','= 20','= 30','= 40','= 50','= 60','= 70','= 80','= 90','= 100'}),'Orientation','horizontal');
title(lgd2,'Edge variability-based threshold: Percentile')
xticks([]);
xlim([0.5 80.5]);
ytickangle(90)

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

Alabel = annotation('textbox',[0, 0.8637, 0.0160, 0.0280],'String','A','EdgeColor','none','FontSize',32);

Blabel = annotation('textbox',[0, 0.5881, 0.0160, 0.0280],'String','B','EdgeColor','none','FontSize',32);

Clabel = annotation('textbox',[0, 0.3125, 0.0160, 0.0280],'String','C','EdgeColor','none','FontSize',32);
