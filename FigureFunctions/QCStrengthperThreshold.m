function QCStrengthperThreshold(PSig_con,PSig_var,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS)

% INPUTS:
%
% PSig_con = A matrix of the propotion of nodes with a significant QCStrength
% correlation for a pipeline (rows) under a given edge based consistency 
% threshold (columns)
%
% PSig_var = A matrix of the propotion of nodes with a significant QCStrength
% correlation for a pipeline (rows) under a given edge based variability 
% threshold (columns)
%
% PROCESSING_MATRIX = A matrix indicating which processing choices were
% used
%
% PROCESSING_MATRIX_LABELS = A cell wehere each cell indicates what each
% row of PROCESSING_MATRIX is

n = size(PSig_con,1)+1;

prop_sig_thr = PSig_con(:,2:11);
prop_sig_thr2 = PSig_var(:,3:12);

[X,Y] = meshgrid(1:size(prop_sig_thr,1),1:size(prop_sig_thr,2));

X = X';
Y = Y';

figure('units','pixels','outerposition',[0 0 2560 1440])

subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.05], [0.1 0.01], [0.1 0.01]); 
subplot_tight(88,1,12:38)

for i = 1:10
idx = find(Y==i);
scatter(X(idx),prop_sig_thr(idx),50,Y(idx),'filled');
hold on
end
colormap(flipud(make_cmap('red',10)))
ylim([0 1.1])
box on
lgd1 = legend({'= 5%','= 10%','= 20%','= 30%','= 40%','= 50%','= 60%','= 70%','= 80%','= 90%'},'Orientation','horizontal');
title(lgd1,'Edge consistency-based threshold: Percentage')
xlim([0.5 n-.5])
ytickangle(90)
xticks([]);
ylabel({'Proportion of edges with a','significant QC-Strength correlation'})
set(gca,'FontSize',16);

subplot_tight(88,1,39:65)
for i = 10:-1:1
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
xlim([0.5 n-.5]);
ytickangle(90)

ylabel({'Proportion of edges with a','significant QC-Strength correlation'})

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


