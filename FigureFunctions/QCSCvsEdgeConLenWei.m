function QCSCvsEdgeConLenWei(QCSC_DATA,LENGTH_DATA,CONSISTENCY_DATA,WEIGHT_DATA,PipelineINDS)

% INPUTS:
%
% QCSC_DATA = A cell array where each cell contains a vector of each edges QCSC correlation p value
%
% LENGTH_DATA = A cell array where each cell contains a vector of each
% edges length
%
% CONSISTENCY_DATA = A cell array where each cell contains a vector of each edges consistency
%
% WEIGHT_DATA = A cell array where each cell contains a vector of each
% edges weight
%
% PipelineINDS = A 1*5 array indicating which pipelines to display

cmap1 = flipud(make_cmap('red',100));

figure('units','pixels','outerposition',[0 0 2560 680])
    
for i = 1:length(PipelineINDS)
    idx = PipelineINDS(i);
    
ax1 = subplot(2,5,i);
scatter(log(WEIGHT_DATA{idx}),QCSC_DATA{idx},10,LENGTH_DATA{idx},'filled')
xlabel('Edge weight (log)')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge lnegth';
set(gca,'FontSize',10);
colormap(ax1,cmap1)
title(['Pipeline ',num2str(idx)],'FontSize',14);

ax2 = subplot(2,5,i+5);
scatter(log(WEIGHT_DATA{idx}),QCSC_DATA{idx},10,CONSISTENCY_DATA{idx},'filled')
xlabel('Edge weight (log)')
ylabel('QC-SC correlation')
c = colorbar;
c.Label.String = 'Edge consistency';
set(gca,'FontSize',10);
colormap(ax2,cmap1)
title(['Pipeline ',num2str(idx)],'FontSize',14);

end

for ROW = 1:2
    for COL = 1:5

PlotLabelLetters = {'A','B','C','D','E';'F','G','H','I','J';'K','L','M','N','O';'P','Q','R','S','T'};
PlotLabel = annotation('textbox',[0.095+(0.163*(COL-1))    0.99-(0.48*(ROW-1))    0.0164    0.0378],'String',PlotLabelLetters{ROW,COL},'EdgeColor','none','FontSize',32);

    end
end

end