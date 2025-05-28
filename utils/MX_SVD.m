% function  [phiOut, rankRare, isRare, residuR] = MXSVD_xw(k,X,row,col)
function [j,target,Ind_set,m,utarget,target_loc,utarget_loc,Eigen_T,Eigen_U]= MX_SVD(k,X,row,col)
% MX-SVD algorithm
% k given signal-rank

% initialize
[p,N] = size(X);
%omega = zeros(p,N);
isRare = zeros(1,N);
count = 1;
Ind_set=[];
% first iteration
     [u,s,v] = svd(X,'econ');
     phi(:,1:k) = u(:,1:k);
     aa=inv(phi'*phi);
     XX = X - phi*phi'*X;
     for i = 1:N
         residu(i) = norm(XX(:,i)); %X(:,i) - phi*phi'*X(:,i));
     end
     [r(1),ind] = max(residu);
     %r(1) = temp^2;
     phiT(:,:,1) = phi;
     % normalize it
     omega = X(:,ind);
     isRare(ind) = count;
     count = count+1;
     Ind_set=[Ind_set,ind];
% all iteration
h=waitbar(0,'please wait');
for j = 1:k
    waitbar(j/k,h) 
    XX = X - omega*inv(omega'*omega)*omega'*X;
     [u,s,v] = svd(XX,'econ');
     if (k-j)>=1
        %psi(:,1:k-j) = u(:,1:k-j);
        %phi(:,1:k-j) = psi(:,1:k-j);
        phi(:,1:k-j) = u(:,1:k-j);
     end
     phi(:,k-j+1:k) = omega;
     
     XX = X - phi*inv(phi'*phi)*phi'*X;
     for i = 1:N
         residu(i) = norm(XX(:,i));      %X(:,i) - phi*inv(phi'*phi)*phi'*X(:,i));
         if isRare(i) ~= 0
             residu(i) = -inf;
         end
     end
     [r(j+1),ind] = max(residu);
     Ind_set=[Ind_set,ind];
     %r(j+1) = temp^2;
     phiT(:,:,j+1) = phi;
     omega(:,1:j) = omega;
     % normalize it
     omega(:,j+1) = X(:,ind);
     isRare(ind) = count;
     count = count+1;
end
delete(h);
% output
[residuR,ind] = min(r);
phiOut = phiT(:,:,ind);
rankRare = ind-1;
isRare(isRare>rankRare) = 0;
j=rankRare;
m=k-j;
target=omega(:,1:j);
Eigen_T=phiOut(:,1:j);
utarget=omega(:,j+1:k);
Eigen_U=phiOut(:,j+1:k);
figure()
plot(r(1,1:k));
[target_loc,utarget_loc]=target_show(rankRare,k,Ind_set,X,row,col);
%% show location of signal in Ind_set
function [target_loc,utarget_loc]=target_show(j,p,Ind_set,img0,row,col)
Ind_set=Ind_set';
x_set=fix(Ind_set./row);
y_set=Ind_set-x_set*row;
for k=1:p
    if x_set(k,1)==0
        x_set(k,1)=x_set(k,1)+1;
    end
    if y_set(k,1)==0
        y_set(k,1)=row;
    end
end
im=reshape(img0(30,:),[row,col]);
target_loc=[x_set(1:j,1) y_set(1:j,1)];
utarget_loc=[x_set(j+1:p,1) y_set(j+1:p,1)];

% figure; imagesc(im); colormap(gray); hold on
% axis off
% axis equal
% for m=1:size(x_set(1:j,1),1)
%     plot(x_set(m,1),y_set(m,1),'*','color','r');
%     text(x_set(m,1)+2,y_set(m,1),num2str(m),'color','r','FontSize',12);
% end

% figure; imagesc(im); colormap(gray); hold on
% axis off
% axis equal
% for m=j+1:size(x_set(1:p,1),1)
%     plot(x_set(m,1),y_set(m,1),'*','color','r');
%     text(x_set(m,1)+1,y_set(m,1),num2str(m),'color','r','FontSize',12);
% end

% figure; imagesc(im); colormap(gray); hold on
% axis off
% axis equal
% for m=1:size(x_set(1:p,1),1)
%     plot(x_set(m,1),y_set(m,1),'*','color','r');
%     text(x_set(m,1)+1,y_set(m,1),num2str(m),'color','r','FontSize',12);
% end
end
end
