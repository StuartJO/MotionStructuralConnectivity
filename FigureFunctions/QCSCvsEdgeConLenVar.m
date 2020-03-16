function QCSCvsEdgeConLenVar(QCSC_DATA,LENGTH_DATA,CONSISTENCY_DATA,VARIABILITY_DATA,PipelineINDS)

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

figure('units','pixels','outerposition',[0 0 2560 680])
%figure('units','pixels','outerposition',[0 0 2560 1440])
    
for i = 1:length(PipelineINDS)
    idx = PipelineINDS(i);

ax1 = subplot(2,5,i);
scatter(CONSISTENCY_DATA{idx},QCSC_DATA{idx},10,LENGTH_DATA{idx},'filled')
xlabel('Edge consistency')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge length';
set(gca,'FontSize',10);
colormap(ax1,cmap1)
title(['Pipeline ',num2str(idx)],'FontSize',14);

ax2 = subplot(2,5,i+5);
scatter(-log(VARIABILITY_DATA{idx}),QCSC_DATA{idx},10,LENGTH_DATA{idx},'filled')
xlabel('Edge weight variability (-log)')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge length';
set(gca,'FontSize',10);
colormap(ax2,cmap1)
title(['Pipeline ',num2str(idx)],'FontSize',14);

end

for ROW = 1:2
    for COL = 1:5

PlotLabelLetters = {'A','B','C','D','E';'F','G','H','I','J';'K','L','M','N','O';'P','Q','R','S','T'};
annotation('textbox',[0.095+(0.163*(COL-1))    0.99-(0.48*(ROW-1))    0.0164    0.0378],'String',PlotLabelLetters{ROW,COL},'EdgeColor','none','FontSize',32);

    end
end

end