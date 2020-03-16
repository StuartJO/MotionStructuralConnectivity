function W = threshold_arbmeasure(W, M, p)
%THRESHOLD_ARBMEASURE     Threshold edges ranked by given arbitrary measure
%
%   W_thr = threshold_arbmeasure(W, M, p);
%
%   This function "thresholds" the connectivity matrix by preserving a
%   proportion p (0<p<1) of the edges with the largest values of the given
%   measure M. All other weights, and all weights on the main diagonal
%   (self-self connections) are set to 0. 
%
%   Inputs: W,      weighted or binary connectivity matrix
%           M,      matrix of values of the measure by which to rank each
%                   node, preserving only edges with the largest values
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights removed)
%
%   Output: W_thr,  thresholded connectivity matrix
%
%   James Roberts, QIMR Berghofer, 2014
%   Based directly on BCT's threshold_proportional.m by:
%     Mika Rubinov, U Cambridge,
%     Roan LaPlante, Martinos Center, MGH
%    (BCT: brain-connectivity-toolbox.net)

n=size(W,1);                                %number of nodes
W(1:n+1:end)=0;                             %clear diagonal

if isequal(W,W.');                          %if symmetric matrix
    W=triu(W);                              %ensure symmetry is preserved
    ud=2;                                   %halve number of removed links
else
    ud=1;
end

ind=find(W);                                %find all links
E=sortrows([ind M(ind)], -2);               %sort by M
en=round((n^2-n)*p/ud);                     %number of links to be preserved

W(E(en+1:end,1))=0;                         %apply threshold

if ud==2                                    %if symmetric matrix
    W=W+W.';                                %reconstruct symmetry
end