function PlotBrainNodalProperty(data,plotName,plotLabel,ColorLim,Colormap,PlotColorBar,ColorbarLim,ignore0)

% This function projects some data onto the cortical surface of fsaverage
%
% INPUTS:
%
% data = a vector of some value for each node
%
% plotName = the title of the plot
%
% plotLabel = the label of the plot (e.g. '1', '2', 'A', 'B' etc etc)
%
% ColorLim = the range of the colormap
%
% Colormap = the colormap to use. Set as 1, 2 or 3 for a blue/red, red, or
% viridisplus colormap respective. Otherwise input an N*3 colormap for use
% (haven't extensively tested this)
%
% PlotColorBar = set to 1 to plot the colorbar (note that plotName becomes
% the titl of the colormap)
%
% ColorbarLim = the range of values to display in the colorbar
%
% ignore0 = set to 1 to color any nodal region with a value of 0 as the
% background colour of the brain ([.5 .5 .5])

runplotLabel = 1;

if nargin < 3
    runplotLabel = 0;
end

if nargin < 4
    ColorLim = [min(data) max(data)];
end

if nargin < 5
   Colormap = 1; 
end

if nargin < 6
    PlotColorBar = 0;
end

if nargin < 7
    ColorbarLim = ColorLim;
end

if nargin < 8
   ignore0 = 0;
end

cmax = max(ColorLim);
cmin = min(ColorLim);

if size(Colormap,2) == 3
    
    cmap = [.5 .5 .5; Colormap];
else
    
if Colormap == 1
    cmap = [.5 .5 .5; make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
elseif Colormap == 2
    cmap = [.5 .5 .5; flipud(make_cmap('orangered',100,30,0))];
elseif Colormap == 3
    cmap = [ .5 .5 .5; viridisplus(100)];
end



end

% This should make it so when the colormaps are combined, the vertices with
% no ROI information should always be coloured grey but everything else
% will be coloured according to the colormap

cmap_length = size(cmap,1)-1;

crange = linspace(cmin,cmax,cmap_length);

grey_val = cmin-diff(crange(1:2));

% Extract the vertex and face information

load('fsaverage_surface_data.mat')
addpath ./cbrewer

lhsurface.vertices = lh_verts;
lhsurface.faces = lh_faces;

rhsurface.vertices = rh_verts;
rhsurface.faces = rh_faces;

% Rotate the data matrix if needed

    if size(data,1) > 1
        data = data';
    end
    
% Remap the values of toeach node onto the respective vertices that make up
% the ROI
    
if length(data) == 220
    cdatal = changem(lh_rand200,[grey_val data(1:100)],0:100);
    cdatar = changem(rh_rand200,[grey_val data(111:210)],0:100);
    
elseif length(data) == 82
    cdatal = changem(lh_aparc,[grey_val data(1:34)],0:34);
    cdatar = changem(rh_aparc,[grey_val data(42:75)],0:34);
    
elseif length(data) == 380
    cdatal = changem(lh_HCPMMP1,[grey_val data(1:180)],0:180);
    cdatar = changem(rh_HCPMMP1,[grey_val data(191:370)],0:180);   
    
end

% If ignoring any ROIS with a value of 0, set their value to grey_val, so
% they will be plotted in grey instead

if ignore0
    cdatal(cdatal == 0) = grey_val;
    cdatar(cdatar == 0) = grey_val;
end

figure('Position',[0 0 1035 791])

subplot_tight = @(m,n,p) subtightplot(m,n,p,[0.005 0.005], [0.01 0.01], [0.01 0.01]); 

subplot_tight(2,2,1)
p = patch(lhsurface);
set(p,'FaceVertexCData',cdatal,'EdgeColor','none','FaceColor','flat');
view([-90 0])
camlight;material dull
caxis([grey_val cmax])
axis off

subplot_tight(2,2,2)
p = patch(rhsurface);
set(p,'FaceVertexCData',cdatar,'EdgeColor','none','FaceColor','flat');
view([90 0])
camlight;material dull
caxis([grey_val cmax])
axis off

subplot_tight(2,2,3)
p = patch(lhsurface);
set(p,'FaceVertexCData',cdatal,'EdgeColor','none','FaceColor','flat');
view([90 0])
camlight;material dull
caxis([grey_val cmax])
axis off

subplot_tight(2,2,4)
p = patch(rhsurface);
set(p,'FaceVertexCData',cdatar,'EdgeColor','none','FaceColor','flat');
view([-90 0])
camlight;material dull
caxis([grey_val cmax])
axis off

% Add the colormap to all plots

colormap(cmap)

% Add either a name to the enter of the plots or a colorbar

if PlotColorBar == 0
    annotation('textbox',[0.3543    0.5228    0.4392    0.0337],'String',plotName,'EdgeColor','none','FontSize',48);
elseif PlotColorBar == 1
    c = colorbar('Location','north');
c.Position = [0.2818    0.4621    0.4392    0.0337];
c.Label.String = plotName;
c.FontSize = 16;
set(c, 'xlim', ColorbarLim)
end

if runplotLabel

annotation('textbox',[0.0117    0.9039    0.0588    0.0885],'String',plotLabel,'EdgeColor','none','FontSize',48);

end
