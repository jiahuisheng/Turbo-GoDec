clear all
close all
warning off
clc    
addpath('utils')
addpath('results')
save_folder = 'results';
       %% Parameter setting
      for data_ind = 1:4 

       %% Read data
          tic
          fprintf('Per-processing VD and NWHFC...');
          [HIM,GT,name] = read_data(data_ind);
          [no_lines,no_rows,no_bands] = size(HIM); 
          img0 = (ToVector(HIM))';% L*N
          fprintf('The Data is');  disp(char(name)); 
       % 

        %% Detector
           method = {'RX-AD','LSMAD','PTA','Auto-AD','RX-AD-BS','LSMAD-BS','PTA-BS','Auto-AD-BS'};



            %% comprehensive compare
           if data_ind == 1 % urban-1
              load urban-1_RXAD_result.mat;load urban-1_LSMAD_result.mat;
              load urban-1_PTA_result.mat;load urban-1_AutoAD_result.mat;
              load urban-1_bs_RXAD_result.mat;load urban-1_bs_LSMAD_result.mat;
              load urban-1_bs_PTA_result.mat;load urban-1_bs_AutoAD_result.mat;
              load urban-1_RGAE_result.mat;load urban-1_bs_RGAE_result.mat;
              AutoAD_result = reshape(AutoAD_result, [no_lines*no_rows, 1]);
              AutoAD_BS_result = reshape(AutoAD_BS_result, [no_lines*no_rows, 1]);
              detection_cub = [RXAD_result(:,1) LSMAD_result(:,1) PTA_result(:) AutoAD_result(:) RXAD_BS_result(:,1) LSMAD_BS_result(:,1) PTA_BS_result(:) AutoAD_BS_result(:)];
           elseif data_ind == 2 % abu-beach-4
              load abu-beach-4_RXAD_result.mat;load abu-beach-4_LSMAD_result.mat;
              load abu-beach-4_PTA_result.mat;load abu-beach-4_AutoAD_result.mat;
              load abu-beach-4_bs_RXAD_result.mat;load abu-beach-4_bs_LSMAD_result.mat;
              load abu-beach-4_bs_PTA_result.mat;load abu-beach-4_bs_AutoAD_result.mat;
              AutoAD_result = reshape(AutoAD_result, [no_lines*no_rows, 1]);
              AutoAD_BS_result = reshape(AutoAD_BS_result, [no_lines*no_rows, 1]);
              detection_cub = [RXAD_result(:,1) LSMAD_result(:,1) PTA_result(:) AutoAD_result(:) RXAD_BS_result(:,1) LSMAD_BS_result(:,1) PTA_BS_result(:) AutoAD_BS_result(:)];
           elseif data_ind == 3 % Hyperion              
              load Hyperion_RXAD_result.mat;load Hyperion_LSMAD_result.mat;
              load Hyperion_PTA_result.mat;load Hyperion_AutoAD_result.mat;
              load Hyperion_bs_RXAD_result.mat;load Hyperion_bs_LSMAD_result.mat;
              load Hyperion_bs_PTA_result.mat;load Hyperion_bs_AutoAD_result.mat;
              AutoAD_result = reshape(AutoAD_result, [no_lines*no_rows, 1]);
              AutoAD_BS_result = reshape(AutoAD_BS_result, [no_lines*no_rows, 1]);
              detection_cub = [RXAD_result(:,1) LSMAD_result(:,1) PTA_result(:) AutoAD_result(:) RXAD_BS_result(:,1) LSMAD_BS_result(:,1) PTA_BS_result(:) AutoAD_BS_result(:)];
           elseif data_ind == 4 % abu-airport-4
              load abu-airport-4_RXAD_result.mat;load abu-airport-4_LSMAD_result.mat;
              load abu-airport-4_PTA_result.mat;load abu-airport-4_AutoAD_result.mat;
              load abu-airport-4_bs_RXAD_result.mat;load abu-airport-4_bs_LSMAD_result.mat;
              load abu-airport-4_bs_PTA_result.mat;load abu-airport-4_bs_AutoAD_result.mat;
              AutoAD_result = reshape(AutoAD_result, [no_lines*no_rows, 1]);
              AutoAD_BS_result = reshape(AutoAD_BS_result, [no_lines*no_rows, 1]);
              detection_cub = [RXAD_result(:,1) LSMAD_result(:,1) PTA_result(:) AutoAD_result(:) RXAD_BS_result(:,1) LSMAD_BS_result(:,1) PTA_BS_result(:) AutoAD_BS_result(:)];
            end               
 

           %% plot-all
            [~] = plot_image(detection_cub,GT,name,method);
            [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SBPR,TDBS,ODP] = Cal_3DROC(detection_cub,GT);
            Criteria = [AUC_pfpd,AUC_tpd,AUC_tpf,TD,BS,SBPR,TDBS,ODP];
            % AUC_AD = round(Criteria,4);
            % excel_name = ['AUC_', num2str(data_ind), '.xlsx'];
            % writematrix(AUC_AD, excel_name);


            % Calculate AUC 2022 version
            AUC_ADP = AUC_tpd;
            AUC_BDP = 1-AUC_tpf;
            AUC_JAD = AUC_pfpd+AUC_ADP;
            AUC_JBS = AUC_pfpd+AUC_BDP;
            AUC_ADBS = AUC_ADP+AUC_BDP;
            AUC_SBPR = AUC_tpd./AUC_tpf;
            AUC_OADP = AUC_pfpd+AUC_ADP+AUC_BDP;
            AUC_AD_new = [AUC_pfpd AUC_tpf AUC_ADP AUC_BDP AUC_JAD AUC_JBS AUC_ADBS AUC_SBPR AUC_OADP];  
            AUC_AD = [AUC_pfpd AUC_ADP AUC_BDP AUC_JAD AUC_JBS AUC_ADBS AUC_SBPR AUC_OADP];
            AUC_AD = round(AUC_AD,4);
            AUC_AD_new = roundn(AUC_AD_new,-4) ; 

            excel_name = ['AUC_new_', num2str(data_ind), '.xlsx'];
            writematrix(AUC_AD, excel_name);
            dlmwrite('Criteria2.txt',char(datetime('today')),'-append','delimiter','\t','precision',5,'newline','pc');
            dlmwrite('Criteria2.txt',char(name),'-append','delimiter','\t','precision',5,'newline','pc');
            dlmwrite('Criteria2.txt',Criteria,'-append','delimiter','\t','precision','%.4f','newline','pc');

            dlmwrite('Criteria_new.txt',char(datetime('today')),'-append','delimiter','\t','precision',5,'newline','pc');
            dlmwrite('Criteria_new.txt',char(name),'-append','delimiter','\t','precision',5,'newline','pc');
            dlmwrite('Criteria_new.txt',AUC_AD_new,'-append','delimiter','\t','precision',5,'newline','pc');
       end
