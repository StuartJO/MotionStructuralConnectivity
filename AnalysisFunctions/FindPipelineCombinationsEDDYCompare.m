function [ORDERED_INDS,ORDERED_MATRIX,LABELS] = FindPipelineCombinationsEDDYCompare(COMBINATIONS,TYPE,ORDER,EXCLUDE)

% This function will pull out the index or indices of a desired pipeline

% INPUTS:
% Required
% COMBINATIONS = a matrix of the combinations of processing choices for
% EDDY1.5 pipelines
%
% TYPE = a 1*7 vector indicating which pipeline/s with a particular 
% parameter you want to extract the index/indices of. The position of each
% value indicates the processing step and the value at that position
% indicates which specific step to extract. A value of zero indicates to
% not extract pipelines based on that processing step. See below for an 
% explanation of how to use
%
% Optional
% ORDER = a 1*7 vector specifying the order in which to arrange the indices 
% (if extracting a single pipeline this has no effect). The first value 
% indicates by which processing step to order pipelines by, the second 
% value indicates which processing step to further order pipelines after
% the first ordering, the third indicates which to use for ordering after
% the second etc etc. 
%
% EXCLUDE = a binary value of 1 or 0. If set to 1, any processing step for 
% which only a single value is being extracted will be removed from 
% ORDERED_MATRIX and LABELS
%
% OUTPUTS:
% ORDERED_INDS = the index/indicies of the desired pipelines, ordering
% according to what is specified by ORDER
%
% ORDERED_MATRIX = a 7*N matrix, where N is the number of pipelines
% extracted. Each column corresponds to a pipeline and each row is a
% processing step. Each row correponds corresponds to a processing step
% (and the ordering of these rows is determined by ORDER), and the value
% indicates the specific algorithm/method being used.
%
% LABELS = a cell indicating what the values in each row of ORDERED_MATRIX
% correspond to. This is set up so if you do imagesc(ORDERED_MATRIX) and
% yticklabels(LABELS) it will format it as per the paper
%
% HOW TO USE:
% TYPE = [MotionCorr TractAlgor SptlCons SeedAlgor TractReWei EdgeWei Parcellation]
%
% MotionCorr = 1, EDDY1 
% MotionCorr = 2, EDDY2 
%
% TractAlgor = 1, FACT 
% TractAlgor = 2, iFOD2 
%
% SptlCons = 1, ACT 
% SptlCons = 2, GWM 
%
% SeedAlgor = 1, dynamic 
% SeedAlgor = 2, WM 
% SeedAlgor = 3, GMWMI 
%
% Filtering = 1, None
% Filtering = 2, SIFT2 
%
% EdgeWei = 1, SSW
% EdgeWei = 2, FA 
%
% Parcellation = 1, 82 Node 
% Parcellation = 2, 220 Node 
% Parcellation = 3, 380 Node 
%
% So TYPE = [1 1 2 3 1 1 3] would extract the index of the pipeline that 
% used EDDY1, FACT, ACT, GMWMI seeding, no filtering, SSW edge weighting
% and the 380 node parcellation.
%
% TYPE = [1 0 0 0 0 0 2] would extract the indicies of all pipelines using
% the 220 node parcellation, TYPE = [1 0 0 0 0 0 2] would extract all
% pipelines using EDDY1 and the 220 node parcellation etc.
%
% ORDER = [7 1 6 2 4 3 5] extracts indices in the order of first being
% sorted by 'Parcellation', then 'MotionCorr', then 'EdgeWei', then
% 'TractAlgor', then 'SptlCons', then 'SeedAlgor', and then by 'TractReWei'.  
% In other words, pipelines are primarily ordered by the type of
% parcellation they used, then within that ordering they are order by the
% type of motion corretion used etc etc
%
% if TYPE = [0 0 1 0 0 0 0], ORDER = [7 1 6 2 4 3 5] and EXCLUDED = 0, then 
% ORDERED_MATRIX(1,:) and LABELS{1} will correspond to the values of
% 'EdgeWei', ORDERED_MATRIX(2,:) and LABELS{2} will correspond to the
% values of 'SptlCons' etc etc. If EXCLUDED = 1, then ORDERED_MATRIX(1,:),  
% LABELS{1} will correspond to the values of 'EdgeWei', ORDERED_MATRIX(2,:),
% LABELS{2} will correspond to the values of 'SeedAlgor' etc etc

if nargin < 2
    % Default ordering
    ORDER = [7 1 6 2 4 3 5];
end

if nargin < 3
    EXCLUDE = 0;
end

EXTRACT_INDS = find(TYPE);

N_TYPES_TO_EXTRACT = length(EXTRACT_INDS);

if N_TYPES_TO_EXTRACT == 0
    
    INDS = 1:size(COMBINATIONS,1);

elseif N_TYPES_TO_EXTRACT == 1

  INDS = find(COMBINATIONS(:,EXTRACT_INDS) == TYPE(EXTRACT_INDS));
  
else

    for i = 1:N_TYPES_TO_EXTRACT
        EXTRACT_IND = EXTRACT_INDS(i);
INDStemp = find(COMBINATIONS(:,EXTRACT_IND) == TYPE(EXTRACT_IND));

if i > 1
   INDS = union(INDS,INDStemp);
else
   INDS = INDStemp;
end

    end
    
end

% I flip the order so that the first processing step pipelines are ordered
% by is the last position in the labels naming convection

ORDER_FLIP = flip(ORDER);

% Because COMBINATION contains a row with values of 1 and 3 instead of 1
% and 3, we turn that 3 into a 2

ORDERED_MATRIX = COMBINATIONS(INDS,:)';

ORDERED_MATRIX = ORDERED_MATRIX(ORDER_FLIP,:);

orderIND = cell(1,7);

for i = 1:7
   [~,orderIND{i}] = sort(ORDERED_MATRIX(i,:));
    ORDERED_MATRIX = ORDERED_MATRIX(:,orderIND{i}); 
end
order = orderIND{6}(orderIND{7});
for i = 5:-1:1
   order =  orderIND{i}(order);
end

ORDERED_INDS = INDS(order);

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;

LABELS = cell(7,1);

LABELS{(ORDER_FLIP==1)} = [sprintf('MotionCorr:{\\color[rgb]{%f,%f,%f}EDDY1}/',Color1),sprintf('{\\color[rgb]{%f,%f,%f}EDDY1.5}/',Color3),sprintf('\\color[rgb]{%f,%f,%f}EDDY2',Color2)];
LABELS{(ORDER_FLIP==6)} = [sprintf('EdgeWei:{\\color[rgb]{%f,%f,%f}SSW}/',Color1),sprintf('\\color[rgb]{%f,%f,%f}FA',Color2)];
LABELS{(ORDER_FLIP==2)} = [sprintf('TractAlgor:{\\color[rgb]{%f,%f,%f}FACT}/',Color1),sprintf('\\color[rgb]{%f,%f,%f}iFOD2',Color2)];
LABELS{(ORDER_FLIP==4)} = [sprintf('SeedAlgor:{\\color[rgb]{%f,%f,%f}dynamic}/',Color1),sprintf('{\\color[rgb]{%f,%f,%f}WM}/',Color2),sprintf('{\\color[rgb]{%f,%f,%f}GMWMI}',Color3)];
LABELS{(ORDER_FLIP==3)} = [sprintf('SptlCons:{\\color[rgb]{%f,%f,%f}ACT}/',Color1),sprintf('\\color[rgb]{%f,%f,%f}GWM',Color2)];
LABELS{(ORDER_FLIP==5)} = [sprintf('TractReWei:{\\color[rgb]{%f,%f,%f}None}/',Color1),sprintf('\\color[rgb]{%f,%f,%f}SIFT2',Color2)];
LABELS{(ORDER_FLIP==7)} = [sprintf('Parcellation:{\\color[rgb]{%f,%f,%f}82 Nodes}/',Color1),sprintf('{\\color[rgb]{%f,%f,%f}220 Nodes}/',Color2),sprintf('{\\color[rgb]{%f,%f,%f}380 Nodes}',Color3)];

if EXCLUDE
    
    EXCLUDE_INDS = TYPE~=0;

    LABELS(EXCLUDE_INDS) = [];

    ORDERED_MATRIX(EXCLUDE_INDS,:) = [];

end
