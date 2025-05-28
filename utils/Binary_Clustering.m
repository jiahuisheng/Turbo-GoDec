function Binary_Clustering(HIM,GT)
% Plot spectral 
img = (ToVector(HIM))';
[m,n,~] = size(HIM);
GT_vector = reshape(GT,[1,m*n]); 
Ind_t = find(GT_vector==1);
Ind_b = find(GT_vector==0);

% target
Target = img(:,Ind_t');
figure()
plot(Target)
mean_c1 = mean(Target');
figure()
plot(mean_c1)

% BKG
BKG=img(:,Ind_b');
% figure()
% plot(BKG)
mean_c2=mean(BKG');
% figure()
% plot(mean_c2)

t_b = [mean_c1' mean_c2'];
figure()
plot(t_b);
legend('Target','BKG');
fprintf('The Euclidean between target and BKG is'); 
pdist2(mean_c1,mean_c2, 'Euclidean') %'mahalanobis'   'cosine'  'Euclidean'  'chebychev' 
fprintf('The cosine between target and BKG is'); 
pdist2(mean_c1,mean_c2, 'cosine') %'mahalanobis'   'cosine'  'Euclidean'  'chebychev'  
