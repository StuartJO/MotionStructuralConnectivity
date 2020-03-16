function [adjGrp adjThr adjMask adjProp] = connectomeGroupThreshold(mats, grpThr, denSD) 

% function [adjGrp adjThr adjMask adjProp] = connectomeGrpThreshold(mats, grpThr, denSD) 

% This function will apply a threshold to each individual connectivity
% matrix based on the proportion of subjects containing a non-zero value in
% each edge. Specifically, for each edge, it determines the proportion of
% subjects with a non-zero value, and retains only those edges whose
% proportion exceeds a threshold value. This is loosely based on the
% approach used in van den Heuvel and Sporns (2011) J Neurosci.
%
% -------
% INPUTS:
% -------
%
% mats  - an array containing K cells each with an N*N matrix, where K =
%         number of subjects and N = number of nodes. Could optionally be a
%         3d matrix of dimensions N*N*K.
%
% grpThr - Group-level threshold. This is the proportion of subjects who
%          must have a non-zero value for a given edge to be retained. Must
%          be between zero and one.
%
% denSD  - Some subjects may have a low connection density and bias the
%          results. This input allows one to exclude those subjects when
%          considering the proportion threshold, grpThr. This value should
%          represent the number of SDs above/below the mean for a person to
%          be excluded from consideration. e.g., denSD = 2 means that any
%          subject with a connectiont density +/- 2SD from mean will be
%          excluded. Note: the subject is only excluded from computing the
%          proportion of subjects with a non-zero edge value. The data from
%          the subject will still be retained and the threshold applied to
%          this subject in the outputs. To use the whole sample, set denSD
%          = 0. This is the default.
%
% -------
% OUTPUTS:
% -------
%
% adjGrp   - a group averaged connectome. The weights are the average
%          across subjects, considering non-zero values only (as per van den Heuvel
%          and Sporns).
%
% adjThr   - an N*N*K matrix of thresholded connectomes for each person.
%
% adjMask  - the binary mask of edges used to threshold connectomes.
%
% adjProp  - Weights of this matrix indicate the the proportion of subjects 
%           containing a non-zero edge (excluding subjects excluded through
%           denSD)
%
% Alex Fornito, Monash University, Jun 2016.
%==========================================================================

%--------------------------------------------------------------------------
% Preliminaries
%--------------------------------------------------------------------------

% set default denSD
if nargin<3
    denSD = 0;
end

% combine data into 3d matrix if input is a cell array
if iscell(mats)
    for i = 1:length(mats)  
        adj(:,:,i) = mats{i};
    end  
else
    adj = mats;
end

% number of subjects
Nsub = size(adj,3);

% number of nodes
Nnodes = size(adj, 1);

%--------------------------------------------------------------------------
% Select edges to retain
%--------------------------------------------------------------------------

% if a matrix has NaN values, raplace them with zero; 
%if any(isnan(adj))~=0
adj(isnan(adj))=0;
%end

% binarize matrices and compute proportion of subjects with each edge

adjBin = adj;
adjBin(adjBin~=0) = 1;

% if option chosen to exclude subjects with low density from calculating
% proportions
if denSD ~= 0 
    
    % compute densities
    for j = 1:size(adjBin, 3)
        den(j) = density_und(adjBin(:,:,j));
    end
    
    % indicies of subjects to exclude
    inds1 = find(den>mean(den)+(denSD*std(den)));
    inds2 = find(den<mean(den)-(denSD*std(den)));
    
    % remove subjects from consideration
    adjBin(:,:,[inds1 inds2]) = [];
    
    Nsub = Nsub - length([inds1 inds2]); 
    
end

% Proportion of subjects with a non-zero edge
adjProp = sum(adjBin,3)./Nsub;

% generate binary 'network mask' of edges to retain
adjMask = adjProp;
adjMask(adjMask<grpThr) = 0;
adjMask(adjMask~=0) = 1;

% print density of group connectome
%fprintf('Density of group connectome is %d \n',density_und(adjMask));

% threshold individual matrices
for k = 1:size(adj,3)
    adjThr(:,:,k) = adj(:,:,k).*adjMask;
end

%--------------------------------------------------------------------------
% Generate group mean matrix (average weights over all non-zero edges)
%--------------------------------------------------------------------------

% initalize output
adjGrp = zeros(Nnodes, Nnodes);

% get average of nonzero vals for each edge in upper triangle
for l = 1:Nnodes
    
    for j = l+1:Nnodes
        
        vals = nonzeros(adjThr(l,j,:));
        
        if ~isempty(vals)
            adjGrp(l, j) = nanmean(vals);
        end
        
    end
end

% symmetrize network
adjGrp = adjGrp+adjGrp';
                