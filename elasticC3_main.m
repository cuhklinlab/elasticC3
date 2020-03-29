function [clu_X_mat_re,clu_Y_mat_re,clu_Z1_mat_re,clu_Z2_mat_re,obj1_mat_re,obj2_mat_re,jsdxy_mat_re,jsdz_mat_re,clu_X_bes,clu_Y_bes] = elasticC3_main(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer, niter_inter, alpha1, alpha2, iter,dist)

%%main function for elasticC3 algorithm
%%initialize the input for niter_outer times
clu_X_mat = zeros(niter_outer,size(X,1));
clu_Y_mat = zeros(niter_outer,size(Y,1));
clu_Z_mat = zeros(niter_outer,size(X,2));
for l = 1:niter_outer
clu_X_mat(l,:) = randsample(nclu_X, size(X,1),true);
clu_Y_mat(l,:) = randsample(nclu_Y, size(Y,1),true);
clu_Z_mat(l,:) = randsample(nclu_Z, size(X,2),true);
end

%%generate empty lists in order to save results obtained by each outer
% For each matrix inside the lists, the first row is the initialization of
% clustering;
clu_X_mat_re = zeros(niter_inter+1,size(X,1),niter_outer);% The last row in each matrix is the result after lable switch
clu_Y_mat_re = zeros(niter_inter+1,size(Y,1),niter_outer);
clu_Z1_mat_re = zeros(niter_inter+1,size(X,2),niter_outer);
clu_Z2_mat_re = zeros(niter_inter+1,size(X,2),niter_outer);
obj1_mat_re = zeros(niter_inter,niter_outer);
obj2_mat_re = zeros(niter_inter,niter_outer);
jsdxy_mat_re = zeros(niter_inter,niter_outer);
jsdz_mat_re = zeros(niter_inter,niter_outer);

%%outer iteration
title1 = sprintf('elasticC3 for NO.%d instance is in progress',iter);
title2 = sprintf('remaining time of elasticC3 for NO.%d instance =',iter);

    
h = waitbar(0,title1);
s = clock;
for i = 1:niter_outer
    %begin process
    [clu_X_Trace, clu_Y_Trace, clu_Z1_Trace,clu_Z2_Trace, obj_val_Trace1, obj_val_Trace2, JSD_Trace_XY, JSD_Trace_Z] = elasticC3(X, Y, clu_X_mat(i,:), clu_Y_mat(i,:), clu_Z_mat(i,:), nclu_X, nclu_Y, nclu_Z, niter_inter, alpha1, alpha2,dist);
    clu_X_mat_re(:,:,i) = clu_X_Trace;
    clu_Y_mat_re(:,:,i) = clu_Y_Trace;
    clu_Z1_mat_re(:,:,i) = clu_Z1_Trace;
    obj1_mat_re(:,i) = obj_val_Trace1;
    jsdxy_mat_re(:,i) = JSD_Trace_XY;
    clu_Z2_mat_re(:,:,i) = clu_Z2_Trace;
    obj2_mat_re(:,i) = obj_val_Trace2;
    jsdz_mat_re(:,i) = JSD_Trace_Z;
    %end process
    %begin estimate remaining time
    if i == 1
        is = etime(clock,s);
        esttime = is*niter_outer;
    end
    h = waitbar(i/niter_outer,h,...
    [title2,num2str(esttime-etime(clock,s),'%4.1f'),'sec' ]);
    %end estimate remaining time
end
close(h);

% choose the best result and return
obj1 = obj1_mat_re(niter_inter,:);
bes1 = find(obj1 == min(obj1), 1 );
clu_X_bes = clu_X_mat_re(niter_inter,:,bes1);
obj2 = obj2_mat_re(niter_inter,:);
bes2 = find(obj2 == min(obj2), 1 );
clu_Y_bes = clu_Y_mat_re(niter_inter,:,bes2);