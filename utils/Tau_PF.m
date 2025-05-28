function [AUC]=Tau_PF(Tao,PF_hydice,w,integral_mode)
%% defination
name = {'RX-AD','LSMAD','PTA','Auto-AD','RX-AD-BS','LSMAD-BS','PTA-BS','Auto-AD-BS'};
line_type = {'k','g','m',[0.2 0.5 0.7],[0.75,0.75,0],[0 0.5 0],'r',"#FF8800",'b'};

%% load data
% figure()
% if w==1
%    plot(Tao(:,1),PF_hydice(:,1),'-','Color',[line_type{1}]','LineWidth',2);
% else
tao_min_max = 0;
if w>1
    figure()
    for i=1:w
        tao=Tao(:,i);pf=PF_hydice(:,i);
        tao_min = min(tao(tao>0));
        if tao_min > tao_min_max
            tao_min_max = tao_min;
        end
        plot(tao,pf,'-','Color',[line_type{i}]','LineWidth',1.5);
        hold on
    end   

    grid on
    xlim([tao_min_max, 1])
    set(gca,'xscale','log')
    set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
%     title ('P_{F} vs \tau')
    xlabel('\tau','FontAngle', 'italic');ylabel('P_{f}','FontAngle', 'italic');
    if w==1
        legend (name{1});
    else
        legend (name{1},name{2},name{3},name{4},name{5},name{6},name{7},name{8},'Box', 'off'); 
    end
end
%% calculate AUC
for j=1:w    
        switch integral_mode
        case 1
    % case 1--matlab version
        AUC0=trapz(Tao(:,j),PF_hydice(:,j));
        case 2
    % case 2--step function
        tmp_Tao=Tao(:,j);
        diet_tmp_Tao=diff(tmp_Tao);
        tmp_PF=PF_hydice(:,j);
        AUC0=sum(diet_tmp_Tao.*tmp_PF(1:end-1));
        end
%***********************************************
    b=PF_hydice(1,j);
    AUC(j,1)=AUC0/b;
end
    