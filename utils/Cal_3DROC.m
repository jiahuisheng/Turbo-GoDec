function [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP] = Cal_3DROC(detection_cub,Ground_Truth)
[N,w] = size(detection_cub); 
Th_val = 1;
integral_mode = 2;
Th_mode = 2;
detection_map = detection_cub;
%**************************************************** GT mode
% mode 1 (three columns)
[row,col] = size(Ground_Truth);
gt = reshape(Ground_Truth,[1,N]); 
% gt(gt~=0)=1;
% % mode 2 (the third column)
% gt=zeros(m,n);
% for  k=1:5
%      col3=[13 23 33 43 53];
%      [loc_row,loc_col]=find(GT==col3(1,k));
%      gt(loc_row,loc_col)=1;
% end
%      gt=reshape(gt,[1,N]); 
%****************************************************
%****************************************************
pd = zeros(N,w);pf = zeros(N,w);Tau = zeros(N,w);
for k = 1:w
    mat = abs(detection_map(:,k));
    min_ds = min(mat);
    max_ds = max(mat);
    mat = (mat-min_ds)/(max_ds-min_ds);
%****************************************************Tau val
switch Th_val
    case 1
    Tau_tmp = sort(mat);N_Tau=N;
    case 2
    Tau_tmp = [0:0.01:1]';N_Tau=101;
end
Tau(:,k) = Tau_tmp;
% pd=zeros(N_Tau,w);pf=zeros(N_Tau,w);
%*****************************************************
for Tindx = 1:N_Tau
    matb = mat;
    th = Tau(Tindx,k);
   %**************************************************Tau mode
   % VERSION 1
%    for Th_mode=2
       switch Th_mode
           case 1  % 
               matb(mat>th) = 1;                
               matb(mat<=th) = 0;
           case 2
               matb(mat>=th) = 1;                
               matb(mat<th) = 0;
           case 3
               matb(mat>=th) = 1;                
               matb(mat<th) = 0;
       end
   
      [~,y] = find(gt>0);
      D = length(find(matb(y,1)==1));
      F = sum(sum(matb))-D;     
      pd(Tindx,k) = D/sum(gt);
      pf(Tindx,k) = F/(N-sum(gt)); 
%    end
end
end
%***********************************3D ROC
  [AUC_pfpd] = PF_PD(pf,pd,w,integral_mode);
  [AUC_tpd] = Tau_PD(Tau,pd,w,integral_mode);
  [AUC_tpf] = Tau_PF(Tau,pf,w,integral_mode);            
  rdD_ROC(Tau,pf,pd,w);
  
  ODP = AUC_pfpd+AUC_tpd-AUC_tpf;
  TD = AUC_pfpd+AUC_tpd;
  BS = AUC_pfpd-AUC_tpf;
  TDBS = AUC_tpd-AUC_tpf;
  SNPR = AUC_tpd./AUC_tpf;
end      
