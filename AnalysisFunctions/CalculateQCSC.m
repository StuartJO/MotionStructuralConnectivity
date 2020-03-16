function [DATA,STRDATA] = CalculateQCSC(motion,edgeconthr,calc_node_STRs,include_zero_edges,save_edge_properties,save_output_location)

% This function calculates the correlation between edge weight and some
% measure of motion across participants.
%
% INPUTS:
%
% motion = a number from 1 to 7, or the name of the motion measure. Choose
% from 'RMS_abs' (1),'RMS_rel' (2),'meanABSrms_b0' (3),'meanABSrms_b3000' 
% (4),'meanRELrms_b0' (5),'meanRELrms_b3000' (6), 'tsnr_b3000' (7)
%
% edgeconthr = edge consistency-based threshold to use (see note below)
%
% calc_node_STRs = set to 1 to calculate node strengths and degrees
%
% include_zero_edges = set to 1 to include edges that have no connections
% in the QCSC correlation calculation
% 
% save_edge_properties = set to 1 to save edge properties in a seperate
% file
%
% save_output_location = set the location of where output is saved
% (defaults to current directory)
%
% OUTPUTS:
%
% DATA = a structure of QCSC related information and edge properties. See 
% below for explanation of all fields in this structure. Each contains a 
% cell/vector of some property for each pipeline. 
%
% STRDATA = a structure of all data related to node degree and strength. 
% See below for explanation of all fields in this structure. Each contains
% a cell/vector of some property for each pipeline. 
% 
%
% Note if calc_node_STRs is set to 1, a group adjacency matrix mask 
% will be created under using a threshold set by 'edgeconthr', which is  
% then applied to each subjects adjacency matrix. These thresholded 
% matrices are further thresholded at a number of different levels to get 
% strength and degree values for the nodes. These two levels of 
% thresholding will have no effect on edge consistency based thresholding, 
% but will effect edge variability based thresholding

if nargin < 6
    save_output_location = [];
end

load('MOTION_DATA.mat','motion_data','MOTIONNAMES','subjects')
load('COMBINATIONS_MATRIX.mat')

if isstring(motion)

switch motion
    
    case 'RMS_abs'
        m = 1;
    case 'RMS_rel'
        m = 2;
    case 'meanABSrms_b0'
        m = 3;
    case 'meanABSrms_b3000'
        m = 4;
    case 'meanRELrms_b0'
        m = 5;
    case 'meanRELrms_b3000'
        m = 6;
    case 'tsnr_b3000'
        m = 7;
end

else
    m = motion;
end

NPipes = size(COMBINATIONS,1);

% Note that DATA1 and DATA2 are combined into DATA, they are just kept
% seperate so they can be easily saved

% A cell array where each cell contains a matrix of each edges QCSC correlation value
DATA1.EdgeMatQCSC = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges QCSC correlation p value
DATA1.EdgeMatQCSC_PVALS = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges consistency
DATA2.EdgeMatConsistency = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges weight variability
DATA2.EdgeMatWeightVariability = cell(NPipes,1);

% A vector of the mean QCSC correlation for each pipeline
DATA1.mean_QCSC = zeros(NPipes,1);

% A vector of the median QCSC correlation for each pipeline
DATA1.median_QCSC = zeros(NPipes,1);

% A cell array where each cell contains a vector of each edges QCSC correlation value
DATA1.QCSC = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges QCSC correlation p value
DATA1.QCSC_PVALS = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges consistency
DATA2.EdgeConsistency = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges weight variability
DATA2.EdgeWeightVariability = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges length
DATA2.EdgeLength = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges length
DATA2.EdgeMatLength = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges weight covariance
DATA2.EdgeCovariance = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges weight variance
DATA2.EdgeWeightVariance = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges weight covariance
DATA2.EdgeMatCovariance = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges weight variance
DATA2.EdgeMatWeightVariance = cell(NPipes,1);

% A cell array where each cell contains a vector of each edges weight
DATA2.EdgeWeight = cell(NPipes,1);

% A cell array where each cell contains a matrix of each edges weight
DATA2.EdgeMatWeight = cell(NPipes,1);

% A vector of the mean total strength (across participants) for each pipeline
DATA2.total_str_mean = zeros(NPipes,1);

% A vector of the standard deviation of total strength (across participants) for each pipeline
DATA2.total_str_sd = zeros(NPipes,1);

% A vector of the mean density (across participants) for each pipeline
DATA2.density_mean = zeros(NPipes,1);

% A vector of the standard deviation of density (across participants) for each pipeline
DATA2.density_sd = zeros(NPipes,1);

% A vector of the mean of mean edge weight (across participants) for each pipeline
DATA2.MeanEdgeWeight_mean = zeros(NPipes,1);

% A vector of the mean of standard deviation of edge weight (across participants) for each pipeline
DATA2.MeanEdgeWeight_sd = zeros(NPipes,1);

% A vector of the standard deviation of mean edge weight (across participants) for each pipeline
DATA2.SdEdgeWeight_mean = zeros(NPipes,1);

% A vector of the standard deviation of standard deviation of edge weight (across participants) for each pipeline
DATA2.SdEdgeWeight_sd = zeros(NPipes,1);

% Thresholds to check when calculating strength
threshs = [0 .05 .1:.1:1];

Nthr = length(threshs);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based consistency threshold (columns)
STRDATA.STRcon = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based consistency threshold (columns)
STRDATA.DEGcon = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns)
STRDATA.STRvar = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns)
STRDATA.DEGvar = cell(NPipes,Nthr);

% A cell matrix where each cell is the QCStrength correlation for each node 
% for a pipeline (rows) under a given edge based consistency threshold (columns)
STRDATA.QCSTR_con = cell(NPipes,Nthr);

% A cell matrix where each cell is the QCStrength correlation p value for each node 
% for a pipeline (rows) under a given edge based consistency threshold (columns)
STRDATA.QCSTR_pval_con = cell(NPipes,Nthr);

% A cell matrix where each cell is the QCStrength correlation for each node 
% for a pipeline (rows) under a given edge based variability threshold (columns)
STRDATA.QCSTR_var = cell(NPipes,Nthr);

% A cell matrix where each cell is the QCStrength correlation p value for each node 
% for a pipeline (rows) under a given edge based variability threshold (columns)
STRDATA.QCSTR_pval_var = cell(NPipes,Nthr);

% A matrix of the propotion of nodes with a significant QCStrength
% correlation for a pipeline (rows) under a given edge based consistency 
% threshold (columns)
STRDATA.PropNodeStr_con_sig = zeros(NPipes,Nthr);

% A matrix of the propotion of nodes with a significant QCStrength
% correlation for a pipeline (rows) under a given edge based variability 
% threshold (columns)
STRDATA.PropNodeStr_var_sig = zeros(NPipes,Nthr);

STRDATA.threshs = threshs; 

% Rarely, when calculating FA MRtrix can encounter an issue
% which results in a NaN value. This seems to be the result of the FA mask
% being slightly smaller than the brain mask, thus streamlines can exist
% just outside of FA values. When assigning FA values to streamlines and
% subsequently edge weights, this will result in a NaN value. We just set
% this value to 0 thereby excluding the edge for that subject. To be very
% conservative, we can also exclude an edge if for any subject a NaN is
% recorded. Doing so does not fundamentally change the results at all.

IgnoreEdgesWithNaN = 0;

for ITER = 1:size(COMBINATIONS,1)
    
    %% Extract data and compute DATA1.QCSC
tic
clear adjs lens con
    load(['Pipeline_',num2str(ITER),'.mat'],'adjs','lens','con')
            
            if m == 1 || m == 2
                movement_data = motion_data{m}(:,COMBINATIONS(ITER,1));
            else
                movement_data = motion_data{m};
            end
    
    
            Nsubs = length(adjs);
            NNodes = length(con);
                
             groupMask = con;
             groupMask(groupMask<edgeconthr) = 0;
             groupMask(groupMask~=0) = 1;
             
                        groupMaskTriu = triu(groupMask,1);
                        
                        [indEdge] = find(groupMaskTriu(:)>0);
                        
                        NEdges = length(indEdge);
                        
                        subjEdges = zeros(Nsubs, NEdges); 
                        subjEdgesLength = zeros(Nsubs, NEdges); 
                        
                        Ws = zeros(NNodes,NNodes,Nsubs);
                        den = zeros(Nsubs,1);
                        total_STR = zeros(Nsubs,1);
                        EdgeWeight_mean = zeros(Nsubs,1);
                        EdgeWeight_sd = zeros(Nsubs,1);
                        
                        for subj=1:Nsubs
                            
                            adjVect = adjs{subj}(:);
                            lengthVect = lens{subj}(:);
                            
                            edge_weights = adjVect(indEdge);
                                                        
                            subjEdges(subj,:) = edge_weights;                            
                            
                            subjEdgesLength(subj,:) = lengthVect(indEdge);
                            W = adjs{subj}.*groupMask;
                                                        
                            Ws(:,:,subj) = W;
                            
                            Nnodes = size(W,1);
                            Nedges = nnz(W(~isnan(W)));
                            den(subj,1) = Nedges/((Nnodes^2-Nnodes)/2);
                            
                            total_STR(subj,1) = nansum(edge_weights);
                            EdgeWeight_mean(subj,1) = nanmean(edge_weights(edge_weights>0));
                            EdgeWeight_sd(subj,1) = nanstd(edge_weights(edge_weights>0));                          
                            
                        end

                        if IgnoreEdgesWithNaN == 1
                            NaN_edges = isnan(sum(subjEdges));
                            subjEdges(:,NaN_edges) = [];
                            subjEdgesLength(:,NaN_edges) = [];
                            
                            NEdges = size(subjEdges,2);
                        end
                        
                        if ~include_zero_edges
                            
                            subjEdges(subjEdges==0) = nan;
                            
                        end
                        
                        DATA2.total_str_mean(ITER,1) = mean(total_STR);
                        DATA2.total_str_sd(ITER,1) = std(total_STR);
                        
                        DATA2.density_mean(ITER,1) = mean(den);
                        DATA2.density_sd(ITER,1) = std(den);
                                                
                        DATA2.MeanEdgeWeight_mean(ITER,1) = mean(EdgeWeight_mean);
                        DATA2.MeanEdgeWeight_sd(ITER,1) = std(EdgeWeight_mean);
                        
                        DATA2.SdEdgeWeight_mean(ITER,1) = mean(EdgeWeight_sd);
                        DATA2.SdEdgeWeight_sd(ITER,1) = std(EdgeWeight_sd);
                        
                        
			if calc_node_STRs   

                        threshs = [0 .05 .1:.1:1];
                        
                        for thr = 1:length(threshs)
                            Str = zeros(Nsubs, NNodes); 
                            Deg = zeros(Nsubs, NNodes); 
                            Str2 = zeros(Nsubs, NNodes); 
                            Deg2 = zeros(Nsubs, NNodes); 
                            [~, ~, AdjMaskCon] = connectomeGroupThreshold(adjs,threshs(thr)); 
                            AdjMaskVar = threshold_consistency(Ws, threshs(thr));
                            for subj=1:Nsubs
                            Str(subj,:) = nansum(adjs{subj}.*AdjMaskCon);
                            Deg(subj,:) = nansum(double((adjs{subj}.*AdjMaskCon) > 0));
                            Str2(subj,:) = nansum(adjs{subj}.*AdjMaskVar);
                            Deg2(subj,:) = nansum(double((adjs{subj}.*AdjMaskVar) > 0));
                            end
                            STRDATA.STRcon{ITER,thr} = Str;
                            STRDATA.DEGcon{ITER,thr} = Deg;
                            STRDATA.STRvar{ITER,thr} = Str2;
                            STRDATA.DEGvar{ITER,thr} = Deg2;
                            
                            [STRDATA.QCSTR_con{ITER,thr},STRDATA.QCSTR_pval_con{ITER,thr}] = corr(Str,movement_data,'type','Spearman');
        
                            [STRDATA.QCSTR_var{ITER,thr},STRDATA.QCSTR_pval_var{ITER,thr}] = corr(Str2,movement_data,'type','Spearman');
     
                            STRDATA.PropNodeStr_con_sig(ITER,thr) = sum(STRDATA.QCSTR_pval_con{ITER,thr} < .05)./length(STRDATA.QCSTR_pval_con{ITER,thr});
                            STRDATA.PropNodeStr_var_sig(ITER,thr) = sum(STRDATA.QCSTR_pval_var{ITER,thr} < .05)./length(STRDATA.QCSTR_pval_var{ITER,thr});
                            STRDATA.threshs = threshs; 
     
                        end
            else
                  STRDATA = [];      
			end
                        
                                                
                        qcsc = zeros(1,NEdges);
                        qcsc_pval = zeros(1,NEdges);
                        Edge_Con = zeros(1,NEdges);
                        Edge_Var = zeros(1,NEdges);
                        Edge_Mat_Rho = zeros(NNodes);
                        Edge_Mat_Pval = zeros(NNodes);
                        Edge_Con_Mat = zeros(NNodes);
                        Edge_WeiVariability_Mat = zeros(NNodes);
                        Edge_Length = zeros(1,NEdges);
                        Edge_Length_Mat = zeros(NNodes);
                        Edge_Cov = zeros(1,NEdges);
                        Edge_Weight_variance = zeros(1,NEdges);
                        Edge_Cov_Mat = zeros(NNodes);
                        Edge_Weight_variance_mat = zeros(NNodes);
                        Edge_Weight = zeros(1,NEdges);
                        Edge_Weight_Mat = zeros(NNodes);
                        
                        for edge = 1:NEdges
                            
                            Edges2Corr = subjEdges(:,edge);
                            
                            movement2corr = movement_data;
                            
                            lengths2use = subjEdgesLength(:,edge);
                            
                            nsubs = length(subjEdges(:,edge));
                                                        
                            movement2corr(isnan(subjEdges(:,edge))) = [];
                            
                            Edges2Corr(isnan(subjEdges(:,edge))) = [];
                            
                            lengths2use(isnan(subjEdges(:,edge))) = [];
                                               
                            Edgecovvar = cov(Edges2Corr,movement2corr);
                            
                            Edge_Cov(edge) = Edgecovvar(1,2);
                            
                            Edge_Weight_variance(edge) = Edgecovvar(1,1);
                            
                            Edge_Cov_Mat(indEdge(edge)) = Edge_Cov(edge);
                            
                            Edge_Weight_variance_mat(indEdge(edge)) = Edgecovvar(1,1);

                            Edge_Con(edge) = sum(Edges2Corr~=0)/nsubs;
                            
                            Edge_Var(edge) = std(Edges2Corr) / mean(Edges2Corr);
                            
                            [qcsc(edge),qcsc_pval(edge)] = corr(movement2corr,Edges2Corr,'Type','Spearman');
                            Edge_Mat_Rho(indEdge(edge)) = qcsc(edge);
                            Edge_Mat_Pval(indEdge(edge)) = qcsc_pval(edge);
                            Edge_Con_Mat(indEdge(edge)) = Edge_Con(edge);
                            Edge_WeiVariability_Mat(indEdge(edge)) = Edge_Var(edge);
                            
                            % If we are including subjects with an edge
                            % weight of 0 in the analysis, this can
                            % serverely distort the calculation of edge
                            % length. So I only calculate the mean length
                            % based on subjects who actually have an edge
                            % present
                            Edge_Length(edge) = mean(lengths2use(lengths2use>0));
                                                        
                            Edge_Length_Mat(indEdge(edge)) = Edge_Length(edge);
                            Edge_Weight(edge) = mean(Edges2Corr);
                            Edge_Weight_Mat(indEdge(edge)) = Edge_Weight(edge);
                            
                        end
                        
                        %% Store output in cells
                        
                        DATA1.EdgeMatQCSC{ITER} = Edge_Mat_Rho + Edge_Mat_Rho';
                        DATA1.EdgeMatQCSC_PVALS{ITER} = Edge_Mat_Pval + Edge_Mat_Pval';
                        DATA2.EdgeMatConsistency{ITER} = Edge_Con_Mat + Edge_Con_Mat';
                        DATA2.EdgeMatWeightVariability{ITER} = Edge_WeiVariability_Mat + Edge_WeiVariability_Mat';
                        DATA1.mean_QCSC(ITER,1) = nanmean(qcsc);
                        DATA1.median_QCSC(ITER,1) = nanmedian(qcsc);
                        DATA1.QCSC{ITER} = qcsc;
                        DATA1.QCSC_PVALS{ITER} = qcsc_pval;
                        DATA2.EdgeConsistency{ITER} = Edge_Con;
                        DATA2.EdgeWeightVariability{ITER} = Edge_Var;
                        DATA2.EdgeLength{ITER} = Edge_Length;
                        DATA2.EdgeMatLength{ITER} = Edge_Length_Mat + Edge_Length_Mat';
                        DATA2.EdgeCovariance{ITER} = Edge_Cov;
                        DATA2.EdgeWeightVariance{ITER} = Edge_Weight_variance;
                        DATA2.EdgeMatCovariance{ITER} = Edge_Cov_Mat + Edge_Cov_Mat';
                        DATA2.EdgeMatWeightVariance{ITER} = Edge_Weight_variance_mat + Edge_Weight_variance_mat;
                        DATA2.EdgeWeight{ITER} = Edge_Weight;
                        DATA2.EdgeMatWeight{ITER} = Edge_Weight_Mat + Edge_Weight_Mat';   
                        
                        timetaken = toc;
                        fprintf('Completed %d/240 in %.4f seconds\n',ITER,timetaken)
            
end

DATA1.COMBINATIONS = COMBINATIONS;
DATA1.movement_data = motion_data{m};

DATA2.COMBINATIONS = COMBINATIONS;

if calc_node_STRs && ~isempty(save_output_location)
          
   save([save_output_location,'/Node_degree_strength_thr_',num2str(edgeconthr),'_inc0Edges_',num2str(include_zero_edges),'.mat'],'-struct','STRDATA','-v7.3')

end

if save_edge_properties && ~isempty(save_output_location)
    
    savename = [save_output_location,'/Pipelines_EdgeProperties_thr_',num2str(edgeconthr),'_inc0Edges_',num2str(include_zero_edges),'.mat'];
    disp(['Saved as ',savename])
    save(savename,'-v7.3','-struct','DATA2')

end


if ~isempty(save_output_location)
    
    savename = [save_output_location,'/Pipelines_QCSC_thr_',num2str(edgeconthr),'_',MOTIONNAMES{m},'_inc0Edges_',num2str(include_zero_edges),'.mat'];
    disp(['Saved as ',savename])
    save(savename,'-v7.3','-struct','DATA1')

end

DATA = MergeStructs(DATA1,DATA2);

function [merged_struct] = MergeStructs(struct_a,struct_b)
%%if one of the structres is empty do not merge
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end
%%insert struct a
merged_struct=struct_a;
%%insert struct b
size_a=length(merged_struct);
for j=1:length(struct_b)
    f = fieldnames(struct_b);
    for i = 1:length(f)
        merged_struct(size_a+j).(f{i}) = struct_b(j).(f{i});
    end
end