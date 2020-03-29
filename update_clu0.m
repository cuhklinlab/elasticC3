function [clu_X_temp0, obj] = update_clu0(p_XZ, clu_X, clu_Z, nclu_X, nclu_Z)
clu_X_temp0 = clu_X;
for i = 1:size(p_XZ,1)
    KL_X = zeros(nclu_X,1);
    for k = 1:nclu_X
        clu_X_temp0(i) = k;
        p_XZ_tilde_temp0 = cal_clu_prob(p_XZ, clu_X_temp0, clu_Z, nclu_X, nclu_Z);
        p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ, p_XZ_tilde_temp0, clu_X_temp0, clu_Z);
        p_post_zx = p_XZ(i,:)/sum(p_XZ(i,:));
        p_tilde_zx_tilde = p_tilde_XZ_temp0(i,:)/sum(p_XZ(i,:));
        ind = p_XZ(i,:)~=0;
        KL_X(k) = sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind))));
    end
    clu_X_temp0(i) = find(KL_X == min(KL_X), 1 );
end

p_XZ_tilde_temp0 = cal_clu_prob(p_XZ, clu_X_temp0, clu_Z, nclu_X, nclu_Z);
p_tilde_XZ_temp0 = cal_coclu_prob(p_XZ, p_XZ_tilde_temp0, clu_X_temp0, clu_Z);
obj = 0;
for j = 1:size(p_XZ,1)
    p_post_zx = p_XZ(j,:)/sum(p_XZ(j,:));
    p_tilde_zx_tilde = p_tilde_XZ_temp0(j,:)/sum(p_XZ(j,:));
    ind = p_XZ(j,:)~=0;
    obj = obj + sum(p_post_zx(ind).*log(p_post_zx(ind)./(p_tilde_zx_tilde(ind))));
end

