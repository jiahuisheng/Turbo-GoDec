w = 4;
name = {'HYDICE Urban','Pavia','Hyperion','abu-beach-1'};
line_type = {'b','b',"b",'b'};
% figure()
K = [80, 90, 100, 110, 120];
m = [20, 25, 30, 35, 40];
AUC_K = [0.9981 0.9979 0.9979 0.9980 0.9983;
        0.9990 0.9988 0.9992 0.9775 0.9776;
        0.9989 0.9990 0.9990 0.9991 0.9991;
        0.9981 0.9984 0.9986 0.9987 0.8418];
AUC_m = [0.9983 0.9981 0.9979 0.9973 0.9966;
        0.9778 0.9777 0.9992 0.9990 0.9984;
        0.9990 0.9991 0.9990 0.9990 0.9990;
        0.8419 0.9987 0.9986 0.9983 0.9981];

for i=1:w
    figure()
    tao=K;pf=AUC_K(i,:);
    plot(tao,pf,'-.o','Color',[line_type{i}]','LineWidth',1.5);
    hold on
    grid on
    ylim([0.5, 1])
% set(gca,'xscale','log')
set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
%     title ('P_{F} vs \tau')
xlabel('K');ylabel('AUC_{(D,F)}');
end   

% grid on
% % xlim([tao_min_max, 1])
% % set(gca,'xscale','log')
% set(gca,'FontSize', 16,'FontWeight', 'bold','FontName', 'Times New Roman')
% %     title ('P_{F} vs \tau')
% xlabel('\tau','FontAngle', 'italic');ylabel('AUC_{(D,F)}','FontAngle', 'italic');