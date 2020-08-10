function [DATA,STRDATA] = QCSC(adjs,motion,lens,edgeconthr,calc_node_STRs,include_zero_edges,save_edge_properties,save_output_location)

% This function calculates the correlation between edge weight and some
% measure of motion across participants.
%
% INPUTS:
%
% adjs = a cell array where each cell contains a subjects weighted connectivity matrix
%
% motion = some measure of head motion
%
% lens = a cell array where each cell contains a subjects distance matrix (length of each connection)
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

if nargin < 8
    save_output_location = [];
end

% Note that DATA and DATA are combined into DATA, they are just kept
% seperate so they can be easily saved

% DATA.EdgeMatQCSC: a matrix of each edges QCSC correlation value

% DATA.EdgeMatQCSC_PVALS: a matrix of each edges QCSC correlation p value

% DATA.EdgeMatConsistency: a matrix of each edges consistency

% DATA.EdgeMatWeightVariability: a matrix of each edges weight variability

% DATA.mean_QCSC: the mean QCSC correlation

% DATA.median_QCSC: the median QCSC correlation

% DATA.QCSC: a vector of each edges QCSC correlation value

% DATA.QCSC_PVALS: a vector of each edges QCSC correlation p value

% DATA.EdgeConsistency: a vector of each edges consistency

% DATA.EdgeWeightVariability: a vector of each edges weight variability

% DATA.EdgeLength: a vector of each edges length

% DATA.EdgeMatLength: a matrix of each edges length

% DATA.EdgeCovariance: a vector of each edges weight covariance

% DATA.EdgeWeightVariance: a vector of each edges weight variance

% DATA.EdgeMatCovariance: a matrix of each edges weight covariance

% DATA.EdgeMatWeightVariance: a matrix of each edges weight variance

% DATA.EdgeWeight: a vector of each edges weight

% DATA.EdgeMatWeight: a matrix of each edges weight

% DATA.total_str_mean: the mean total strength (across participants)

% DATA.total_str_sd: standard deviation of total strength (across participants)

% DATA.density_mean: the mean density (across participants)

% DATA.density_sd: the standard deviation of density (across participants) 

% DATA.MeanEdgeWeight_mean: the mean of mean edge weight (across participants)

% DATA.MeanEdgeWeight_sd: the mean of standard deviation of edge weight (across participants)

% DATA.SdEdgeWeight_mean: the standard deviation of mean edge weight (across participants)

% DATA.SdEdgeWeight_sd: the standard deviation of standard deviation of edge weight (across participants)

% Thresholds to check when calculating strength
threshs = [0 .05 .1:.1:1];

Nthr = length(threshs);

% STRDATA.STRcon: the node strength 
% under a given edge based consistency threshold (columns)
STRDATA.STRcon = cell(1,Nthr);

% A cell matrix where each cell is the node degree 
% under a given edge based consistency threshold (columns)
STRDATA.DEGcon = cell(1,Nthr);

% A cell matrix where each cell is the node strength 
% under a given edge based variability threshold (columns)
STRDATA.STRvar = cell(1,Nthr);

% A cell matrix where each cell is the node degree 
% under a given edge based variability threshold (columns)
STRDATA.DEGvar = cell(1,Nthr);

% A cell matrix where each cell is the QCStrength correlation for each node 
%  under a given edge based consistency threshold (columns)
STRDATA.QCSTR_con = cell(1,Nthr);

% A cell matrix where each cell is the QCStrength correlation p value for each node 
%  under a given edge based consistency threshold (columns)
STRDATA.QCSTR_pval_con = cell(1,Nthr);

% A cell matrix where each cell is the QCStrength correlation for each node 
%  under a given edge based variability threshold (columns)
STRDATA.QCSTR_var = cell(1,Nthr);

% A cell matrix where each cell is the QCStrength correlation p value for each node 
% under a given edge based variability threshold (columns)
STRDATA.QCSTR_pval_var = cell(1,Nthr);

% A matrix of the propotion of nodes with a significant QCStrength
% correlation (p < .05) under a given edge based consistency threshold (columns)
STRDATA.PropNodeStr_con_sig = zeros(1,Nthr);

% A matrix of the propotion of nodes with a significant QCStrength
% correlation (p < .05) under a given edge based variability threshold (columns)
STRDATA.PropNodeStr_var_sig = zeros(1,Nthr);

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
   
tic    
    
            Nsubs = length(adjs);
                
                [~, ~,groupMask] = connectomeGroupThreshold(adjs, edgeconthr); 

                NNodes = length(groupMask);
             
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
                        
                        DATA.total_str_mean = mean(total_STR);
                        DATA.total_str_sd = std(total_STR);
                        
                        DATA.density_mean = mean(den);
                        DATA.density_sd = std(den);
                                                
                        DATA.MeanEdgeWeight_mean = mean(EdgeWeight_mean);
                        DATA.MeanEdgeWeight_sd = std(EdgeWeight_mean);
                        
                        DATA.SdEdgeWeight_mean = mean(EdgeWeight_sd);
                        DATA.SdEdgeWeight_sd = std(EdgeWeight_sd);
                        
                        
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
                            STRDATA.STRcon{thr} = Str;
                            STRDATA.DEGcon{thr} = Deg;
                            STRDATA.STRvar{thr} = Str2;
                            STRDATA.DEGvar{thr} = Deg2;
                            
                            [STRDATA.QCSTR_con{thr},STRDATA.QCSTR_pval_con{thr}] = corr(Str,motion,'type','Spearman');
        
                            [STRDATA.QCSTR_var{thr},STRDATA.QCSTR_pval_var{thr}] = corr(Str2,motion,'type','Spearman');
     
                            STRDATA.PropNodeStr_con_sig(thr) = sum(STRDATA.QCSTR_pval_con{thr} < .05)./length(STRDATA.QCSTR_pval_con{thr});
                            STRDATA.PropNodeStr_var_sig(thr) = sum(STRDATA.QCSTR_pval_var{thr} < .05)./length(STRDATA.QCSTR_pval_var{thr});
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
                            
                            movement2corr = motion;
                            
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
                        
                        %% Store output in structures
                        
                        DATA.EdgeMatQCSC = Edge_Mat_Rho + Edge_Mat_Rho';
                        DATA.EdgeMatQCSC_PVALS = Edge_Mat_Pval + Edge_Mat_Pval';
                        DATA.EdgeMatConsistency = Edge_Con_Mat + Edge_Con_Mat';
                        DATA.EdgeMatWeightVariability = Edge_WeiVariability_Mat + Edge_WeiVariability_Mat';
                        DATA.mean_QCSC = nanmean(qcsc);
                        DATA.median_QCSC = nanmedian(qcsc);
                        DATA.QCSC = qcsc;
                        DATA.QCSC_PVALS = qcsc_pval;
                        DATA.EdgeConsistency = Edge_Con;
                        DATA.EdgeWeightVariability = Edge_Var;
                        DATA.EdgeLength = Edge_Length;
                        DATA.EdgeMatLength = Edge_Length_Mat + Edge_Length_Mat';
                        DATA.EdgeCovariance = Edge_Cov;
                        DATA.EdgeWeightVariance = Edge_Weight_variance;
                        DATA.EdgeMatCovariance = Edge_Cov_Mat + Edge_Cov_Mat';
                        DATA.EdgeMatWeightVariance = Edge_Weight_variance_mat + Edge_Weight_variance_mat;
                        DATA.EdgeWeight = Edge_Weight;
                        DATA.EdgeMatWeight = Edge_Weight_Mat + Edge_Weight_Mat';   
                        
                        timetaken = toc;
                        fprintf('Completed in %.4f seconds\n',timetaken)
            

DATA.motion = motion;