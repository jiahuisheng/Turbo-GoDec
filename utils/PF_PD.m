function [AUC]=PF_PD(PF_hydice,PD_hydice,w,integral_mode)
%% defination
% method = {'RX-AD','CEM-AD','ICASC-RXAD','ICASC-CEMAD','CRD','CDA-RXAD','CDA-CEMAD','LSDM-MoG','PTA','Auto-AD','RGAE'};
name = {'RX-AD','LSMAD','PTA','Auto-AD','RX-AD-BS','LSMAD-BS','PTA-BS','Auto-AD-BS'};
% line_type = {'b','g','m',[0.2 0.5 0.7],[0.75,0.75,0],'k','r',[0.8 0.2 0.5],[0 0.5 0]};
line_type = {'k','g','m',[0.2 0.5 0.7],[0.75,0.75,0],[0 0.5 0],'r',"#FF8800",'b'};
% line_type = {'r','b','g','m',[0.2 0.5 0.7],[0.75,0.75,0],'k',[0.8 0.4 0.6],[0 0.5 0],[1 0.5 0],[0.49,0.18,0.56],[0.52 0.36 0.36]};
% [RXAD_result(:,1) CRD_result(:,4) OSP_GoDec_result(:,1) LDSM_result(:) PTA_result(:) RGAE_result(:) ICASC_result(:,1) CDASC_result(:,1) AEIT_result(:,1)];
%% load data

if w>1
    figure()
    for i = 1:w
        pf = PF_hydice(:,i);pd=PD_hydice(:,i);
        plot(pf,pd,'-','Color',[line_type{i}]','LineWidth',1.5);
        hold on
    end   

    grid on
    set(gca,'xscale','log')
    set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
%     title ('P_{D} vs P_{F}')
    xlabel('P_{f}','FontAngle', 'italic');ylabel('P_{d}','FontAngle', 'italic');
    
    if w==1
        legend (name{1});
    else
        legend (name{1},name{2},name{3},name{4},name{5},name{6},name{7},name{8}, 'Location', 'southeast','Box', 'off');  
    end
end



%% calculate AUC
for j=1:w
    switch integral_mode
        case 1
        % case 1--matlab version
           AUC0 = abs(trapz(PF_hydice(:,j),PD_hydice(:,j)));
        case 2
        % case 2--step function
           tmp_PF = PF_hydice(:,j);
           diet_tmp_PF = diff(tmp_PF(end:-1:1));
           tmp_PD = PD_hydice(:,j);
           tmp_PD = tmp_PD(end:-1:1);
           AUC0 = sum(diet_tmp_PF.*tmp_PD(1:end-1));
    end
%***************************************************
        a = PD_hydice(size(PD_hydice,1),j);
        AUC(j,1) = (AUC0-a)/(1-a);        
end
