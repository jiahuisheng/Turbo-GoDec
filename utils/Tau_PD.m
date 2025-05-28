function [AUC]=Tau_PD(Tao,PD_hydice,w,integral_mode)
%% defination
name = {'RX-AD','LSMAD','PTA','Auto-AD','RX-AD-BS','LSMAD-BS','PTA-BS','Auto-AD-BS'};
% line_type = {'b','g','m',[0.2 0.5 0.7],[0.75,0.75,0],'k','r',[0.8 0.2 0.5],[0 0.5 0]};
line_type = {'k','g','m',[0.2 0.5 0.7],[0.75,0.75,0],[0 0.5 0],'r',"#FF8800",'b'};
%% load data
% figure()
% if w==1
%    plot(Tao(:,1),PD_hydice(:,1),'-','Color',[line_type{1}]','LineWidth',2);
% else
if w>1
    figure()
    for i=1:w
        tao=Tao(:,i);pd=PD_hydice(:,i);
        plot(tao,pd,'-','Color',[line_type{i}]','LineWidth',1.5);
        hold on
    end   

    grid on
    xlim([0.0001, 1])
    set(gca,'xscale','log')
    set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
%     title ('P_{D} vs \tau')
    xlabel('\tau','FontAngle', 'italic');ylabel('P_{d}','FontAngle', 'italic');
    if w==1
        legend (name{1});
    else
        legend (name{1},name{2},name{3},name{4},name{5},name{6},name{7},name{8}, 'Location', 'southwest','Box', 'off'); 
    end
end
%% calculate AUC
for j=1:w   
    switch integral_mode
        case 1
    % case 1--matlab version
        AUC0=trapz(Tao(:,j),PD_hydice(:,j));
        case 2
    % case 2--step function
        tmp_Tao=Tao(:,j);
        diet_tmp_Tao=diff(tmp_Tao);
        tmp_PD=PD_hydice(:,j);
        AUC0=sum(diet_tmp_Tao.*tmp_PD(1:end-1));
    end
%***********************************************
    a=PD_hydice(size(PD_hydice,1),j);
    AUC(j,1)=(AUC0-a)/(1-a);
end