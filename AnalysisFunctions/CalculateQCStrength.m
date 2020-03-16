m = 1;
load('Node_degree_strength_thr_0.05_inc0Edges_0.mat')
load('COMBINATIONS_MATRIX.mat')
load('MOTION_DATA.mat')
for ITER = 1:240
   tic     
            if m == 1 || m == 2
                movement_data = motion_data{m}(:,COMBINATIONS(ITER,1));
            else
                movement_data = motion_data{m};
            end
                        threshs = [0 .05 .1:.1:1];
                        
                        for thr = 1:length(threshs)
                            
                            Str = STRcon{ITER,thr};
                            Str2 = STRvar{ITER,thr};
                            
                            [QCSTR_con{ITER,thr},QCSTR_pval_con{ITER,thr}] = corr(Str,movement_data,'type','Spearman');
        
                            [QCSTR_var{ITER,thr},QCSTR_pval_var{ITER,thr}] = corr(Str2,movement_data,'type','Spearman');
     

                            PropNodeStr_con_sig(ITER,thr) = sum(QCSTR_pval_con{ITER,thr} < .05)./length(QCSTR_pval_con{ITER,thr});
                            PropNodeStr_var_sig(ITER,thr) = sum(QCSTR_pval_var{ITER,thr} < .05)./length(QCSTR_pval_var{ITER,thr});
                            
     
                        end
                                timetaken = toc;
                        fprintf('Completed %d/240 in %.4f seconds\n',ITER,timetaken)           
end


save('Node_degree_strength_thr_0.05_inc0Edges_0.mat','DEGcon','STRcon','DEGvar','STRvar','QCSTR_con','QCSTR_pval_con','QCSTR_var','QCSTR_pval_var','PropNodeStr_con_sig','PropNodeStr_var_sig','-v7.3')
