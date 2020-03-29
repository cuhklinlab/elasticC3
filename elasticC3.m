function [clu_X_Trace, clu_Y_Trace, clu_Z1_Trace,clu_Z2_Trace, obj_val_Trace1, obj_val_Trace2, JSD_Trace_XY, JSD_Trace_Z] = elasticC3(X, Y, clu_X, clu_Y, clu_Z, nclu_X, nclu_Y, nclu_Z, niter, alpha1, alpha2,dist)
p_XZ = X/sum(sum(X));
q_YZ = Y/sum(sum(Y));
p_XZ_tilde = cal_clu_prob(p_XZ,clu_X,clu_Z,nclu_X, nclu_Z);

% initialize traces
clu_X_Trace = zeros(niter+1, size(clu_X,2));
clu_Y_Trace = zeros(niter+1, size(clu_Y,2));
clu_Z1_Trace = zeros(niter+1, size(clu_Z,2));
clu_Z2_Trace = zeros(niter+1, size(clu_Z,2));
obj_val_Trace1 = zeros(niter,1);
obj_val_Trace2 = zeros(niter,1);
JSD_Trace_XY = zeros(niter,1);
JSD_Trace_Z = zeros(niter,1);
p_XZ_tilde_Trace = zeros(size(p_XZ_tilde,1),size(p_XZ_tilde,2), niter+1);
clu_X_Trace(1,:) = clu_X;
clu_Y_Trace(1,:) = clu_Y;
clu_Z1_Trace(1,:) = clu_Z;
p_XZ_tilde_Trace(:,:,1) = p_XZ_tilde;

% Step 1   
for iter1 = 1:niter
    [clu_X_Trace(iter1+1,:), obj1] = update_clu0(p_XZ, clu_X_Trace(iter1,:), clu_Z1_Trace(iter1,:), nclu_X, nclu_Z);
    obj_val_Trace1(iter1) = obj1;
    clu_Z1_Trace(iter1+1,:) = update_clu2(p_XZ,clu_X_Trace(iter1+1,:), clu_Z1_Trace(iter1,:),nclu_X,nclu_Z);
    p_XZ_tilde_Trace(:,:,iter1+1) = cal_clu_prob(p_XZ, clu_X_Trace(iter1+1,:), clu_Z1_Trace(iter1,:), nclu_X, nclu_Z);
end


% Step 2
clu_Z2_Trace(1,:) = clu_Z1_Trace(niter+1,:);
for iter2 = 1:niter
    [clu_Y_Trace(iter2+1,:),obj2,jsd_XY] = update_clu_tranf(q_YZ, clu_Y_Trace(iter2,:), clu_Z2_Trace(iter2,:), p_XZ_tilde_Trace(:,:,niter+1), nclu_Y, nclu_Z, alpha1,alpha2,dist);
    [clu_Z2_Trace(iter2+1,:),~,jsd_Z] = update_clu_tranf_Z(q_YZ,clu_Y_Trace(iter2+1,:), clu_Z2_Trace(iter2,:),p_XZ_tilde_Trace(:,:,niter+1),nclu_Y,nclu_Z,alpha2,dist);
    obj_val_Trace2(iter2) = obj2;
    JSD_Trace_XY(iter2) = jsd_XY;
    JSD_Trace_Z(iter2) = jsd_Z;
end

end  