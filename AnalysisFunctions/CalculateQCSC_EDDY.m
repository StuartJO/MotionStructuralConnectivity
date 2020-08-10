function DATA = CalculateQCSC_EDDY(threshold,include_zero_edges,save_output_location)

if nargin < 3
    save_output_location = [];
end

load('MOTION_DATA_EDDY1.5.mat','motion_data')

motion_data_abs = motion_data(:,1);

NPipes = 8;

DATA.EdgeMatQCSC = cell(NPipes,1);
DATA.EdgeMatQCSC_PVALS = cell(NPipes,1);
DATA.EdgeMatConsistency = cell(NPipes,1);
DATA.EdgeMatWeightVariability = cell(NPipes,1);
DATA.mean_QCSC = zeros(NPipes,1);
DATA.median_QCSC = zeros(NPipes,1);
DATA.QCSC = cell(NPipes,1);
DATA.QCSC_PVALS = cell(NPipes,1);
DATA.EdgeConsistency = cell(NPipes,1);
DATA.EdgeWeightVariability = cell(NPipes,1);
DATA.EdgeCovariance = cell(NPipes,1);
DATA.EdgeWeightVariance = cell(NPipes,1);
DATA.EdgeMatCovariance = cell(NPipes,1);
DATA.EdgeMatWeightVariance = cell(NPipes,1);
DATA.EdgeWeight = cell(NPipes,1);
DATA.EdgeMatWeight = cell(NPipes,1);

DATA.total_str_mean = zeros(NPipes,1);
DATA.total_str_sd = zeros(NPipes,1);

DATA.density_mean = zeros(NPipes,1);
DATA.density_sd = zeros(NPipes,1);

DATA.MeanEdgeWeight_mean = zeros(NPipes,1);
DATA.MeanEdgeWeight_sd = zeros(NPipes,1);

DATA.SdEdgeWeight_mean = zeros(NPipes,1);
DATA.SdEdgeWeight_sd = zeros(NPipes,1);

DATA.COMBINATIONS = zeros(NPipes,7);

% Rarely, when calculating FA MRtrix can encounter an issue
% which results in a NaN value. This seems to be the result of the FA mask
% being slightly smaller than the brain mask, thus streamlines can exist
% just outside of FA values.

% When assigning FA values to streamlines and
% subsequently edge weights, this will result in a NaN value. We just set
% this value to 0 thereby excluding the edge for that subject. To be very
% conservative, we can also exclude an edge if for any subject a NaN is
% recorded. Doing so does not fundamentally change the results at all.

IgnoreEdgesWithNaN = 0;

for ITER = 1:8
    
    %% Extract data and compute QCSC
tic
clear adjs con
    load(['Pipeline_',num2str(ITER),'_EDDY1.5.mat'],'adjs','con','COMBINATION')   
    
            Nsubs = length(adjs);
            NNodes = length(con);
    
             groupMask = con;
             groupMask(groupMask<threshold) = 0;
             groupMask(groupMask~=0) = 1;

                        groupMaskTriu = triu(groupMask,1);
                        
                        [indEdge] = find(groupMaskTriu(:)>0);
                        
                        NEdges = length(indEdge);
                        
                        subjEdges = zeros(Nsubs, NEdges);
                        
                        Ws = zeros(NNodes,NNodes,Nsubs);
                        den = zeros(Nsubs,1);
                        total_strength = zeros(Nsubs,1);
                        EdgeWeight_mean = zeros(Nsubs,1);
                        EdgeWeight_sd = zeros(Nsubs,1);
                        
                        for subj=1:Nsubs
                            
                            adjVect = adjs{subj}(:);
                            
                            edge_weights = adjVect(indEdge);
                            
                            subjEdges(subj,:) = edge_weights;      
                            W = adjs{subj}.*groupMaskTriu;
                                                        
                            Ws(:,:,subj) = W;
                            Nnodes = size(W,1);
                            Nedges = nnz(W(~isnan(W)))/2;
                            den(subj,1) = Nedges/((Nnodes^2-Nnodes)/2);
                            
                            total_strength(subj,1) = nansum(edge_weights);
                            EdgeWeight_mean(subj,1) = nanmean(edge_weights(edge_weights>0));
                            EdgeWeight_sd(subj,1) = nanstd(edge_weights(edge_weights>0));                          
                            
                        end
                                              
                        % A later step excluded any participants who have a
                        % NaN value from the analysis. By setting the zeros
                        % to NaN, we exclude participants who do not
                        
                        if IgnoreEdgesWithNaN == 1
                            NaN_edges = isnan(sum(subjEdges));
                            subjEdges(:,NaN_edges) = [];
                            
                            NEdges = size(subjEdges,2);
                        end
                        
                        
                        if ~include_zero_edges
                            
                            subjEdges(subjEdges==0) = nan;
                            
                        end
                        
                        DATA.total_str_mean(ITER,1) = mean(total_strength);
                        DATA.total_str_sd(ITER,1) = std(total_strength);
                        
                        DATA.density_mean(ITER,1) = mean(den);
                        DATA.density_sd(ITER,1) = std(den);
                                                
                        DATA.MeanEdgeWeight_mean(ITER,1) = mean(EdgeWeight_mean);
                        DATA.MeanEdgeWeight_sd(ITER,1) = std(EdgeWeight_mean);
                        
                        DATA.SdEdgeWeight_mean(ITER,1) = mean(EdgeWeight_sd);
                        DATA.SdEdgeWeight_sd(ITER,1) = std(EdgeWeight_sd);
                                         
                        qcsc = zeros(1,NEdges);
                        qcsc_pval = zeros(1,NEdges);
                        Edge_Con = zeros(1,NEdges);
                        Edge_Var = zeros(1,NEdges);
                        Edge_Mat_Rho = zeros(NNodes);
                        Edge_Mat_Pval = zeros(NNodes);
                        Edge_Con_Mat = zeros(NNodes);
                        Edge_WeiVariability_Mat = zeros(NNodes);
                        Edge_Cov = zeros(1,NEdges);
                        Edge_Weight_variance = zeros(1,NEdges);
                        Edge_Cov_Mat = zeros(NNodes);
                        Edge_Weight_variance_mat = zeros(NNodes);
                        Edge_Weight = zeros(1,NEdges);
                        Edge_Weight_Mat = zeros(NNodes);
                        
                        for edge = 1:NEdges
                            
                            Edges2Corr = subjEdges(:,edge);
                            
                            movement2corr = motion_data_abs;
                            
                            nsubs = length(subjEdges(:,edge));
                                                        
                            movement2corr(isnan(subjEdges(:,edge))) = [];
                            
                            Edges2Corr(isnan(subjEdges(:,edge))) = [];
                                                                          
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
                            
                            Edge_Weight(edge) = mean(Edges2Corr);
                            Edge_Weight_Mat(indEdge(edge)) = Edge_Weight(edge);
                            
                        end
                        
                        %% Store output in cells
                        
                        DATA.EdgeMatQCSC{ITER} = Edge_Mat_Rho + Edge_Mat_Rho';
                        DATA.EdgeMatQCSC_PVALS{ITER} = Edge_Mat_Pval + Edge_Mat_Pval';
                        DATA.EdgeMatConsistency{ITER} = Edge_Con_Mat + Edge_Con_Mat';
                        DATA.EdgeMatWeightVariability{ITER} = Edge_WeiVariability_Mat + Edge_WeiVariability_Mat';
                        DATA.mean_QCSC(ITER,1) = nanmean(qcsc);
                        DATA.median_QCSC(ITER,1) = nanmedian(qcsc);
                        DATA.QCSC{ITER} = qcsc;
                        DATA.QCSC_PVALS{ITER} = qcsc_pval;
                        DATA.EdgeConsistency{ITER} = Edge_Con;
                        DATA.EdgeWeightVariability{ITER} = Edge_Var;
                        DATA.EdgeCovariance{ITER} = Edge_Cov;
                        DATA.EdgeWeightVariance{ITER} = Edge_Weight_variance;
                        DATA.EdgeMatCovariance{ITER} = Edge_Cov_Mat + Edge_Cov_Mat';
                        DATA.EdgeMatWeightVariance{ITER} = Edge_Weight_variance_mat + Edge_Weight_variance_mat;
                        DATA.EdgeWeight{ITER} = Edge_Weight;
                        DATA.EdgeMatWeight{ITER} = Edge_Weight_Mat + Edge_Weight_Mat';   
                        DATA.COMBINATIONS(ITER,:) = COMBINATION;
                        
                        timetaken = toc;
                        fprintf('Completed %d/8 in %.4f seconds\n',ITER,timetaken)
            
end

DATA.movement_data = motion_data_abs;

if ~isempty(save_output_location)
    
    savename = [save_output_location,'/Pipelines_EDDY1.5_QCSC_thr_',num2str(threshold),'_ABSall_inc0Edges_',num2str(include_zero_edges),'.mat'];
    disp(['Saved as ',savename])
    save(savename,'-v7.3','-struct','DATA')

end
