function cmap = viridisplus(m)

viridis_data = viridis;

hsv=rgb2hsv(viridis_data);
viridis_resamp=interp1(linspace(0,1,size(viridis_data,1)),hsv,linspace(.1,.98,231));
viridis_resamp=hsv2rgb(viridis_resamp);

warning('off','MATLAB:interp1:UsePCHIP')
YlOrBr_data = cbrewer('seq','YlOrBr',100);

hsv=rgb2hsv(YlOrBr_data);
YlOrBr_resamp=interp1(linspace(0,1,size(YlOrBr_data,1)),hsv,linspace(.25,1,25));
YlOrBr_resamp=hsv2rgb(YlOrBr_resamp);

cm = [viridis_resamp; YlOrBr_resamp];

if nargin < 1
    cmap = cm;
else
    hsv=rgb2hsv(cm);
    cmap=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cmap=hsv2rgb(cmap);
  
end