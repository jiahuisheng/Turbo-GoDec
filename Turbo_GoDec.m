function [L,S,ST,prob_S,RMSE,error]=Turbo_GoDec(X,rank,card,power,row,col)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        GoDec Algotithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%X: nxp data matrix with n samples and p features   
%rank: rank(L)<=rank
%card: card(S)<=card
%power: >=0, power scheme modification, increasing it lead to better
%accuracy and more time cost
%OUTPUTS:
%L:Low-rank part
%S:Sparse part
%RMSE: error
%error: ||X-L-S||/||X||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE:
%Tianyi Zhou and Dacheng Tao, "GoDec: Randomized Lo-rank & Sparse Matrix
%Decomposition in Noisy Case", ICML 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tianyi Zhou, 2011, All rights reserved.

%iteration parameters
iter_max=1e+2;
error_bound=1e-3;
iter=1;
RMSE=[];

%matrix size
[m,n]=size(X);
if m<n 
    X=X'; 
    m = n;
end %must be N*L

%initialization of L and S
L=X;
S=sparse(zeros(size(X)));
% S = full(S);
% rng(0, 'twister')
while true
    %% Update of L
    [U, Sigma, V] = svd(X - S, 'econ');
    L_new = U(:, 1:rank) * Sigma(1:rank, 1:rank) * V(:, 1:rank)';

    %% Update of S
    T=L-L_new+S;
    L=L_new;

    % module B of Turbo-GoDec
    prob_S = Turbo_GoDec_B(T,row,col,10);
    [Temp,idx]=sort(prob_S,'descend');

    S=zeros(size(X));
    new_card = card;
    S(idx(1:new_card),:)=T(idx(1:new_card),:);


    %Error, stopping criteria
    T(idx(1:new_card),:)=0;
    RMSE=[RMSE norm(T(:))];
    if RMSE(end)<error_bound || iter>iter_max
        break;
    else
        L=L+T;
    end
    iter=iter+1;
    
end

LS=L+S;
ST = S+T;
LT = L+T;
error=norm(LS(:)-X(:))/norm(X(:));
if m>n 
    LS=LS';
    L=L';
    S=S';
    ST = ST';
    LT = LT';
end