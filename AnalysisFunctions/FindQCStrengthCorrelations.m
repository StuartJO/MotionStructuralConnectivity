function [C,P] = FindQCStrengthCorrelations(threshold,M)

threshold_values = [0 .05 .1:.1:1];

threshold_ind = find(threshold_values == threshold);

for i = 1:length(STRcon)

    STR = STRcon{i,threshold_ind};

    if M == 1 || M == 2
        
       [C{i},P{i}] = corr(STR,MOTION{M}(subs2use,COMBINATION(i,1)),'type','Spearman');
        
    else
    
        [C{i},P{i}] = corr(STR,MOTION{M}(subs2use,:),'type','Spearman');
    
    end
    
end