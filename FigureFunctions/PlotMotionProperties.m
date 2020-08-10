function PlotMotionProperties

% This function makes some jittered scatters showing the distribution
% ofdifferent motion parameters for the data

subplotind = 1;

load('MOTION_DATA.mat','motion_data','MOTIONNAMES')
figure('Position',[0 0 1280 720])

MotionPlotNames = {'{\itABS}{\it_{all}}','{\itREL}{\it_{all}}','{\itABS}{\it_{b0}}'...
    ,'{\itABS}{\it_{b3000}}','{\itREL}{\it_{b0}}','{\itREL}{\it_{b3000}}','{\itTSNR}'};

for i = 1:7

    if i == 1 || i == 2
        DATA = {motion_data{i}(:,1)};
        subplot(3,3,subplotind)
        JitteredParallelScatter(DATA,1,1,0)
        xticks([])
        ylabel('mm')
        xlabel([MotionPlotNames{i},' [EDDY1]'])
        
        MotionName{subplotind} = [MotionPlotNames{i},' [EDDY1]'];
        meanMotion(subplotind) = mean(DATA{1});    
        stdMotion(subplotind) = std(DATA{1});
        
        
        subplotind = subplotind + 1;
        
        DATA = {motion_data{i}(:,2)};
        %data(:,subplotind) = motion_data{i}(:,2);
        subplot(3,3,subplotind)
        JitteredParallelScatter(DATA,1,1,0)        
        xticks([])
        ylabel('mm')
        xlabel([MotionPlotNames{i},' [EDDY2]'])

        MotionName{subplotind} = [MotionPlotNames{i},' [EDDY2]'];
        meanMotion(subplotind) = mean(DATA{1});    
        stdMotion(subplotind) = std(DATA{1});        

        subplotind = subplotind + 1;
        
    else
        
        DATA = {motion_data{i}};
        %data(:,subplotind) = motion{i};
        subplot(3,3,subplotind)
        JitteredParallelScatter(DATA,1,1,0)
        xticks([])
        xlabel(MotionPlotNames{i})

        
        if i == 7
           ylabel('SNR')
        else
           ylabel('mm') 
        end
        
        MotionName{subplotind} = MOTIONNAMES{i};
        meanMotion(subplotind) = mean(DATA{1});    
        stdMotion(subplotind) = std(DATA{1});
        
        subplotind = subplotind + 1;
        
                
    end
end

PlotLabelLetters = {'A','B','C';'D','E','F';'G','H','I'};

for ROW = 1:3
    for COL = 1:3
        PlotLabel = annotation('textbox',[0.082+(0.28*(COL-1))    0.98-(0.3*(ROW-1))    0.0164    0.0378],'String',PlotLabelLetters{ROW,COL},'EdgeColor','none','FontSize',32);
    end
end