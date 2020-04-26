%% Calculates QCSC values for each edge, and also extracts strength and degree data

% Note this assumes you are in the location where this script is located
PATH = pwd;

addpath(genpath(PATH))

% This calculates the main data (QCSC with ABSall at a .05 edge consistency threshold) reported in the paper

CalculateQCSC(1,.05,1,0,1,'./analysed_data');

% This calculates the main data (QCSC with ABSall) with no threshold applied

CalculateQCSC(1,0,0,1,1,'./analysed_data')

% This calculates QCSC with different motion measures

for i = 2:7
    CalculateQCSC(i,.05,0,0,0,'./analysed_data');
end

% This calculates QCSC for EDDY1.5

CalculateQCSC_EDDY(.05,0,'./analysed_data');
