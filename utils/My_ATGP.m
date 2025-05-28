function [Loc,Sig]=My_ATGP(HIM,M)

bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);

r=reshape(HIM,xx*yy,bnd);
r=r'; 
%=====Find the first point

temp=sum(r.*r);
[a,b]=max(temp);

if (rem(b,xx)==0)
    Loc(1,1)=b/xx;
    Loc(1,2)=xx;
elseif (floor(b/xx)==0)
    Loc(1,1)=1;
    Loc(1,2)=b;
else    
    Loc(1,1)=floor(b/xx)+1;  % y
    Loc(1,2)=b-xx*floor(b/xx);  % x 
end


Sig(:,1)=r(:,b); 
% fprintf('1\n'); 
%==========
for m=2:M
    
    U=Sig;
    P_U_perl=eye(bnd)-U*inv(U'*U)*U';
    y=P_U_perl*r;
    temp=sum(y.*y);
    [a,b]=max(temp);
    if (rem(b,xx)==0)
        Loc(m,1)=b/xx;
        Loc(m,2)=xx;
    elseif (floor(b/xx)==0)
        Loc(m,1)=1;
        Loc(m,2)=b;
    else    
        Loc(m,1)=floor(b/xx)+1;  % y
        Loc(m,2)=b-xx*floor(b/xx);  % x 
    end
    Sig(:,m)=r(:,b);
%     disp(m)
end

% figure; imagesc(HIM(:,:,23)); colormap(gray); hold on
% axis off
% axis equal
% for m=1:size(Loc,1)
%     plot(Loc(m,1),Loc(m,2),'o','color','y');
%     text(Loc(m,1)+2,Loc(m,2),num2str(m),'color','y','FontSize',20);
% end
% 
% figure()
% for tt=1:M
% plot(Sig(:,tt))
% hold on
% end
