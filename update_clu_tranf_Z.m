function [clu_Z_temp, obj, jsd] = update_clu_tranf_Z(p_XZ, clu_X, clu_Z, q_YZ_tilde_temp, nclu_X, nclu_Z, alpha2,dist)
clu_Z_temp = clu_Z; % to store the updated value of the target cluster for each x in X;
for i = 1:size(p_XZ,2)
    KL_Z = zeros(nclu_Z,1);
    for k = 1:nclu_Z
        clu_Z_temp(i) = k;
        p_XZ_tilde_temp = cal_clu_prob(p_XZ, clu_X, clu_Z_temp, nclu_X, nclu_Z);
        p_tilde_XZ_temp = cal_coclu_prob(p_XZ, p_XZ_tilde_temp, clu_X, clu_Z_temp);
        p_post_XZ = p_XZ(:,i)/sum(p_XZ(:,i));
        p_tilde_XZ_tilde = p_tilde_XZ_temp(:,i)/sum(p_XZ(:,i));
        
        ind1 = p_XZ(:,i)~=0;
        temp1 = sum(p_post_XZ(ind1).*log(p_post_XZ(ind1)./(p_tilde_XZ_tilde(ind1))));   
        
         % add the Jesen-Shannon Divergence between Z1 and Z2;
        JSD_Z = dist_cal(p_XZ_tilde_temp,q_YZ_tilde_temp,dist,1);
        KL_Z(k) = temp1*sum(p_XZ(:,i))+ alpha2*JSD_Z;
    end
    clu_Z_temp(i) = find(KL_Z == min(KL_Z), 1 );
end

p_XZ_tilde_temp0 = cal_clu_prob(p_XZ, clu_X, clu_Z_temp, nclu_X, nclu_Z);
p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ, p_XZ_tilde_temp0, clu_X, clu_Z_temp);
obj = 0;
for j = 1:size(p_XZ,1)
    p_post_zx = p_XZ(j,:)/sum(p_XZ(j,:));
    p_tilde_zx_tilde = p_tilde_XZ_temp0(j,:)/sum(p_XZ(j,:));
    ind = p_XZ(j,:)~=0;
    obj = obj + sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind))));
end

jsd = dist_cal(p_XZ_tilde_temp0,q_YZ_tilde_temp,dist,1);


