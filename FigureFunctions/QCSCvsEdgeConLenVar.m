function QCSCvsEdgeConLenVar2(QCSC_DATA,LENGTH_DATA,CONSISTENCY_DATA,VARIABILITY_DATA,PipelineLabels)

% INPUTS:
%
% QCSC_DATA = A cell array where each cell contains a vector of each edges QCSC correlation p value
%
% LENGTH_DATA = A cell array where each cell contains a vector of each
% edges length
%
% CONSISTENCY_DATA = A cell array where each cell contains a vector of each edges consistency
%
% VARIABILITY_DATA = A cell array where each cell contains a vector of each
% edges weight variability
%
% PipelineINDS = A 1*5 array indicating which pipelines to display

cmap1 = flipud(make_cmap('red',100));

figure('units','pixels','outerposition',[0 0 2560 1050])
%figure('units','pixels','outerposition',[0 0 2560 1440])

plot_height = .35;
plot_width = .175;

plot_x_offset = .07;

x_start = .05;
%plot_y_offset = .05;

FONTSIZE = 23;

for i = 1:length(QCSC_DATA)

%ax1 = subplot(2,4,i);

POS_TOP{i} = [x_start+(((plot_x_offset+plot_width)*(i-1))), .6, plot_width, plot_height];
POS_BOTTOM{i} = [x_start+((plot_x_offset+plot_width)*(i-1)), .1, plot_width, plot_height];

ax1 = axes('Position',POS_TOP{i});

scatter(CONSISTENCY_DATA{i},QCSC_DATA{i},10,LENGTH_DATA{i},'filled')
xlabel('Edge consistency')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge length';
set(gca,'FontSize',FONTSIZE);
colormap(ax1,cmap1)
title(PipelineLabels{i},'FontSize',24);
ylim([-1 1])
clims = caxis;
caxis([0 max(clims)])
c.Ticks = 0:50:250;

%ax2 = subplot(2,4,i+4);
ax2 = axes('Position',POS_BOTTOM{i});

scatter(-log(VARIABILITY_DATA{i}),QCSC_DATA{i},10,LENGTH_DATA{i},'filled')
xlabel('Edge weight variability (-log)')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge length';
set(gca,'FontSize',FONTSIZE);
colormap(ax2,cmap1)
%title(PipelineLabels{i},'FontSize',18);
ylim([-1 1])
clims = caxis;
caxis([0 max(clims)])
c.Ticks = 0:50:250;

end

TOPROWLBLS = {'A','B','C','D'};
BOTTOMROWLBLS = {'E','F','G','H'};

LABL_Y_OFFSET = plot_height+.025;

for COL = 1:4
    
    TOP_LBL_POS = [POS_TOP{COL}(1)-.05 POS_TOP{COL}(2)+LABL_Y_OFFSET 0.0164    0.0378];
    annotation('textbox',TOP_LBL_POS,'String',TOPROWLBLS{COL},'EdgeColor','none','FontSize',32);
    
    BOTTOM_LBL_POS = [POS_BOTTOM{COL}(1)-.05 POS_BOTTOM{COL}(2)+LABL_Y_OFFSET 0.0164    0.0378];
    annotation('textbox',BOTTOM_LBL_POS,'String',BOTTOMROWLBLS{COL},'EdgeColor','none','FontSize',32);
    
end

end