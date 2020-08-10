%% Calculates QCSC values for each edge, and also extracts strength and degree data

% Note this assumes you are in the location where this script is located
PATH = pwd;

%addpath(genpath(PATH))

CalculateQCSC(1,.05,1,0,1,'./analysed_data');

CalculateQCSC(1,0,0,1,1,'./analysed_data')

parfor i = [2:9]
    CalculateQCSC(i,.05,0,0,0,'./analysed_data');
end

CalculateQCSC_EDDY(.05,0,'./analysed_data');