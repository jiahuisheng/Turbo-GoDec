function rdD_ROC(Tao,PF_hydice,PD_hydice,w)
%% defination
name = {'RX-AD','LSMAD','PTA','Auto-AD','RX-AD-BS','LSMAD-BS','PTA-BS','Auto-AD-BS'};
% line_type = {'b','g','m',[0.2 0.5 0.7],[0.75,0.75,0],'k','r',[0.8 0.2 0.5],[0 0.5 0]};
line_type = {'k','g','m',[0.2 0.5 0.7],[0.75,0.75,0],[0 0.5 0],'r',"#FF8800",'b'};
%% load data

% if w==1
%         pf=PF_hydice(:,1);tao=Tao(:,1);pd=PD_hydice(:,1);
%         plot3(pf,tao,pd,'-','Color',[line_type{1}]','LineWidth',2);
% else
if w>1
    figure()
    for i=1:w        
        pf=PF_hydice(:,i);tao=Tao(:,i);pd=PD_hydice(:,i);
        plot3(pf,tao,pd,'-','Color',[line_type{i}]','LineWidth',1.5);
        hold on
    end   

    box on
%     set(gca,'xscale','log')
    set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
%     title ('3D ROC')
    xlabel('P_{f}','FontAngle', 'italic');ylabel('\tau','FontAngle', 'italic');zlabel('P_{d}','FontAngle', 'italic');
    if w==1
        legend (name{1});
    else
        legend (name{1},name{2},name{3},name{4},name{5},name{6},name{7},name{8}, 'Location', 'southeast','Box', 'off','FontSize', 10); 
    end
end