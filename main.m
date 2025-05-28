clear all
close all
warning off
clc    
addpath('utils')
addpath('results')
addpath('Dataset')
save_folder = 'results';
method = {'RX-AD','LSMAD','PTA','GTVLRR'};
seed = 2; % 你可以选择任何整数作为种子
rng(seed); % 设置随机数生成器的种子
       %% Parameter setting
         VD_case = 2; % 1--HFC;2--NWHFC
         mj_case = 2; % 1--MX-ATGP; 2--MX-SVD
         th = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];

      for data_ind = 2 % [HU, Pavia, Hyperion, abu-beach-1]     

       %% Read data
          tic
          fprintf('Per-processing VD and NWHFC...');
          [HIM,HIM_norm,GT,name] = read_data(data_ind);
          [no_lines,no_rows,no_bands] = size(HIM); 

          img0_norm = (ToVector(HIM_norm))';% L*N
          img0 = img0_norm;
          fprintf('The Data is');  disp(char(name)); 
       % 
       %% VD calculation
          switch VD_case
             case 1 % HFC
                p = HFC(HIM,0.0001);
                name_VD = ('HFC');

             case 2 % NWHFC
                p = NWHFC(HIM,th(data_ind));
                name_VD = ('NWHFC');
          end
          fprintf('The case of VD is');  disp(char(name_VD)); 

       %% determine j and m       
           switch mj_case
              case 1
                 [j_atgp,target_atgp,Ind_set_atgp,m_atgp,Utarget_atgp,tloc_atgp,utloc_atgp] = MX_ATGP(p,img0,HIM,no_lines,no_rows);                 
                 j = j_atgp;
                 m = p-j_atgp;

              case 2 
                 [j_svd,target_svd,Ind_set_svd,m_svd,Utarget_svd,tloc_svd,utloc_svd,Eigen_T,Eigen_U]= MX_SVD(p,img0,no_lines,no_rows);                 
                 j = j_svd;
                 m = p-j_svd; 
           end
           fprintf('The number of BKG estimated is'); disp(m);
           fprintf('The number of target estimated is'); disp(j);
           t0 = toc;       

           HIM = HIM_norm;

           %% Matrix decomposition
           tic
           fprintf('Turbo-GoDec decomposition processing...');
           % Turbo-GoDec
           [L,S,ST,prob_S,~,~] = Turbo_GoDec(img0',m,j*20,0,no_lines,no_rows); % l*N
           t_TG = toc;


           %% LSMAD_Turbo_GoDec
           tic
           fprintf('Turbo_GoDec method processing...');
           [u_lsmad,~,K_lsmad] = Cal_uRK(L); %  
           
           [LSMAD_result,AUC_LSMAD,AUC_LSMAD_new] = RX_AD(img0,K_lsmad,u_lsmad,no_lines,no_rows,GT,method{2});
           t_TG_1 = toc;


           t_TG_all = t0 + t_TG + t_TG_1;
           disp(t_TG_all)

          
           LSMAD_result_norm = LSMAD_result./max(LSMAD_result(:));
           prob_S_norm = prob_S./max(prob_S(:));
           [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP] = Cal_3DROC(prob_S_norm,GT);
           AUC_AD = [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP]';
           AUC_AD_prob_S = roundn(AUC_AD,-4); 
           S_map = abs(reshape(prob_S_norm,[no_lines,no_rows]));
           % figure(),colormap;imagesc(S_map);axis off;title ('prob S')

           alpha = 0.2;
           s2 = alpha*LSMAD_result_norm+(1-alpha)*prob_S_norm; % urban 0.4 pavia 0.2 hyperion 0.6
           [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP] = Cal_3DROC(s2,GT);
           AUC_AD = [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SNPR,TDBS,ODP]';
           AUC_AD_TG = roundn(AUC_AD,-4); 
           S_map = abs(reshape(s2,[no_lines,no_rows]));
           figure(),colormap;imagesc(S_map);axis off;title ('LSMAD Turbo GoDec+prob S')

           % save (fullfile(save_folder, [char(name) '_TG_result.mat']),'s2');
      end