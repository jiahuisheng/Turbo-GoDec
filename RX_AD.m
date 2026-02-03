function [tmp,AUC_AD,AUC_AD_new] = RX_AD(img0,K0,u0,no_lines,no_rows,GT,algorithm)

N = size(img0,2);
result = (img0-u0*ones(1,N))'*pinv(K0)*(img0-u0*ones(1,N));
 
tmp = diag(result,0);
RXAD_map = abs(reshape(tmp,[no_lines,no_rows]));
% figure(),colormap;imagesc(RXAD_map);axis off;title (algorithm)

% Calculate AUC
[AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP] = Cal_3DROC(abs(tmp),GT);
AUC_AD = [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP]';
AUC_AD = roundn(AUC_AD,-4);  

% Calculate AUC 2022 version
%***********************************3D ROC
AUC_ADP = AUC_tpd;
AUC_BDP = 1-AUC_tpf;
AUC_JAD = AUC_pfpd+AUC_ADP;
AUC_JBS = AUC_pfpd+AUC_BDP;
AUC_ADBS = AUC_ADP+AUC_BDP;
AUC_SNPR = AUC_tpd./AUC_tpf;
AUC_OADP = AUC_pfpd+AUC_ADP+AUC_BDP;
AUC_AD_new = [AUC_pfpd AUC_tpf AUC_ADP AUC_BDP AUC_JAD AUC_JBS AUC_ADBS AUC_SNPR AUC_OADP]';  
AUC_AD_new = roundn(AUC_AD_new,-4);  
