function [prob_S_reshape] = Turbo_GoDec_B(T0,row,col,iteration_LBP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% T0: N by L matrix

% OUTPUTS:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~, L] = size(T0);
% T = sum(abs(T0),2); % T is [0, 1]
% T = T/max(T);
% t = reshape(T, [row, col]);

T = sum(T0,2);
T = T/max(abs(T));
% T = T/std(T);
t_temp = reshape(T, [row, col]);
t = zeros(row+4,col+4);
t(3:row+2, 3:col+2) = t_temp;
N = row+4;
M = col+4;


%% parameter setting
% default 
psi00 = 0.5;
psi01 = 0.3;
psi10 = 0.3;
psi11 = 0.5;

% psi00 = 0.8;
% psi01 = 0.2;
% psi10 = 0.2;
% psi11 = 0.8;

iter_max = 100;


%% initialization
iter = 1;
% N = row;
% M = col;
rl=0.5 * ones(N, M);
rr=0.5 * ones(N, M);
rt=0.5 * ones(N, M);
rb=0.5 * ones(N, M);
% intial_point = 0.2 * ones(row, col);
% rl(3:row+2, 3:col+2) = 0.2 * ones(row, col);
% rr(3:row+2, 3:col+2) = 0.2 * ones(row, col);
% rt(3:row+2, 3:col+2) = 0.2 * ones(row, col);
% rb(3:row+2, 3:col+2) = 0.2 * ones(row, col);

% %% message passing: f->s
% a = sqrt(var1/var0) .* exp( t.^2 ./(2*var1) - t.^2 ./(2*var0) );
% a = 1./ (1 + a);
% pi_in = a;
%% message passing: f->x->g->s
var1 = 0.3;
var2 = 1;
% var1 = 0.8;
% var2 = 1;
tmpVar = var1 + var2;
a = sqrt(tmpVar./var1) .* exp( abs(t).^2 ./tmpVar - abs(t).^2 ./var1 );
a = 1./ (1 + a);
pi_in = a;

%% message passing (Belief Propagation) on 2D Markov chain (RMF)
while(true)
    % step1: update rl, forward message passing on row
    for n = 1:N
        for m = 2:M
            tmp_val = ( psi11*pi_in(n,m-1)*rl(n,m-1)*rt(n,m-1)*rb(n,m-1) + ...
                        psi01*(1-pi_in(n,m-1))*(1-rl(n,m-1))*(1-rt(n,m-1))*(1-rb(n,m-1)) ) ...
                      ./...
                      ( (psi11+psi10)*pi_in(n,m-1)*rl(n,m-1)*rt(n,m-1)*rb(n,m-1) + ...
                        (psi00+psi01)*(1-pi_in(n,m-1))*(1-rl(n,m-1))*(1-rt(n,m-1))*(1-rb(n,m-1)) );
            rl(n,m) = 0.9*rl(n,m)+0.1*tmp_val;
        end
    end

    % step2: update rt, forward message passing on column
    for n = 2:N
        for m = 1:M
            tmp_val = ( psi11*pi_in(n-1,m)*rl(n-1,m)*rt(n-1,m)*rr(n-1,m) + ...
                        psi01*(1-pi_in(n-1,m))*(1-rl(n-1,m))*(1-rt(n-1,m))*(1-rr(n-1,m)) ) ...
                      ./...
                      ( (psi11+psi10)*pi_in(n-1,m)*rl(n-1,m)*rt(n-1,m)*rr(n-1,m) + ...
                        (psi00+psi01)*(1-pi_in(n-1,m))*(1-rl(n-1,m))*(1-rt(n-1,m))*(1-rr(n-1,m)) );
            rt(n,m) = 0.9*rt(n,m)+0.1*tmp_val;
        end
    end

    % step3: update rr, back message passing on row
    for n = N:-1:1
        for m = M-1:-1:1
            tmp_val = ( psi11*pi_in(n,m+1)*rr(n,m+1)*rt(n,m+1)*rb(n,m+1) + ...
                        psi10*(1-pi_in(n,m+1))*(1-rr(n,m+1))*(1-rt(n,m+1))*(1-rb(n,m+1)) ) ...
                      ./...
                      ( (psi11+psi01)*pi_in(n,m+1)*rr(n,m+1)*rt(n,m+1)*rb(n,m+1) + ...
                        (psi00+psi10)*(1-pi_in(n,m+1))*(1-rr(n,m+1))*(1-rt(n,m+1))*(1-rb(n,m+1)) );
            rr(n,m) = 0.9*rr(n,m)+0.1*tmp_val;
        end
    end

    % step4: update rb, back message passing on column
    for n = N-1:-1:1
        for m = M:-1:1
            tmp_val = ( psi11*pi_in(n+1,m)*rl(n+1,m)*rr(n+1,m)*rb(n+1,m) + ...
                        psi10*(1-pi_in(n+1,m))*(1-rl(n+1,m))*(1-rr(n+1,m))*(1-rb(n+1,m)) ) ...
                      ./...
                      ( (psi11+psi01)*pi_in(n+1,m)*rl(n+1,m)*rr(n+1,m)*rb(n+1,m) + ...
                        (psi00+psi10)*(1-pi_in(n+1,m))*(1-rl(n+1,m))*(1-rr(n+1,m))*(1-rb(n+1,m)) );
            rb(n,m) = 0.9*rb(n,m)+0.1*tmp_val;
        end
    end
    
    if iter>iter_max
        break;
    end
    pi_out = (rl.*rr.*rt.*rb)./((1.-rl).*(1.-rr).*(1.-rt).*(1.-rb)+(rl.*rr.*rt.*rb));
    iter = iter+1;
end


%% message passing: s->g->x->f
pi_out = (rl.*rr.*rt.*rb)./((1.-rl).*(1.-rr).*(1.-rt).*(1.-rb)+(rl.*rr.*rt.*rb));


%% cal MAP
% prob_S = (pi_in.*pi_out)./((1.-pi_in).*(1.-pi_out)+(pi_in.*pi_out));
% prob_S_reshape = reshape(prob_S,[row*col,1]);
prob_S_temp = (pi_in.*pi_out)./((1.-pi_in).*(1.-pi_out)+(pi_in.*pi_out));
prob_S = prob_S_temp(3:row+2, 3:col+2);

prob_S_reshape = reshape(prob_S,[row*col,1]);
sda=1;

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% T = sum(abs(T0),2);
% T = sum(T0,2);
% T = T/max(abs(T));
% % T = T/std(T);
% t_temp = reshape(T, [row, col]);
% t = zeros(row+2,col+2);
% t(2:row+1, 2:col+1) = t_temp;
% 
% 
% N = row+2;
% M = col+2;
% % var1 = get_variance_with_percent(T,0.95,'ascend');
% % var2 = get_variance_with_percent(T,0.001,'descend');
% var1 = 0.3;
% var2 = 1;
% 
% iter_max = 10;
% iter = 1;
% 
% pi_in = zeros(N,M);
% 
% rf=zeros(N, M);
% rbr=0.5 * ones(N, M);
% rbl=0.5 * ones(N, M);
% 
% 
% rou000 = 0.997;
% rou001 = 0.003;
% rou010 = 0.8;
% rou011 = 0.2;
% rou100 = 0.8;
% rou101 = 0.2;
% rou110 = 0.3;
% rou111 = 0.7;
% 
% % rf(1,1)=0.071873; % the prob of anomaly
% a = (rou111-rou011-rou101+rou001)^2;
% b = rou011+rou101-2*rou001-1;
% c = rou001;
% rf(1,1)=(-b-sqrt(b^2-4*a*c))/(2*a); % the prob of anomaly
% 
% rou01 = 0.004;
% rou00 = 1-rou01;
% 
% rou10 = (1/rf(1,1)-1)*rou01;
% rou11 = 1-rou10;
% 
% %% message passing: f->x->g->s
% tmpVar = var1 + var2;
% a = sqrt(tmpVar./var1) .* exp( abs(t).^2 ./tmpVar - abs(t).^2 ./var1 );
% a = 1./ (1 + a);
% pi_in = a;
% 
% %% approximate message passing on 2D Markov chain
% while(true)
% % step1: s1,m->h1,m+1->s1,m+1
% for m = 2:M
%     rf(1,m)=(rou01*(1-pi_in(1,m-1)).*(1-rf(1,m-1)).*(1-rbl(1,m-1)) + rou11*pi_in(1,m-1).*rf(1,m-1).*rbl(1,m-1))...
%             ./((1-pi_in(1,m-1)).*(1-rf(1,m-1)).*(1-rbl(1,m-1)) + pi_in(1,m-1).*rf(1,m-1).*rbl(1,m-1));
% end
% 
% % step2: s1,m->h1,m+1->s1,m+1
% for n = 2:N
%     rf(n,1)=(rou01*(1-pi_in(n-1,1)).*(1-rf(n-1,1)).*(1-rbr(n-1,1)) + rou11*pi_in(n-1,1).*rf(n-1,1).*rbr(n-1,1))...
%             ./((1-pi_in(n-1,1)).*(1-rf(n-1,1)).*(1-rbr(n-1,1)) + pi_in(n-1,1).*rf(n-1,1).*rbr(n-1,1));
% end
% 
% % step3: s1,m->h1,m+1->s1,m+1
% for m = M-1:-1:1
%     rbr(1,m)=(rou10*(1-pi_in(1,m+1)).*(1-rbr(1,m+1)).*(1-rbl(1,m+1)) + rou11*pi_in(1,m+1).*rbr(1,m+1).*rbl(1,m+1))...
%             ./((rou00+rou10)*(1-pi_in(1,m+1)).*(1-rbr(1,m+1)).*(1-rbl(1,m+1)) + (rou01+rou11)*pi_in(1,m+1).*rbr(1,m+1).*rbl(1,m+1));
% end
% 
% % step4: s1,m->h1,m+1->s1,m+1
% for n = N-1:-1:1
%     rbl(n,1)=(rou10*(1-pi_in(n+1,1)).*(1-rbr(n+1,1)).*(1-rbl(n+1,1)) + rou11*pi_in(n+1,1).*rbr(n+1,1).*rbl(n+1,1))...
%             ./((rou00+rou10)*(1-pi_in(n+1,1)).*(1-rbr(n+1,1)).*(1-rbl(n+1,1)) + (rou01+rou11)*pi_in(n+1,1).*rbr(n+1,1).*rbl(n+1,1));
% end
% 
% % step5: s1,m->h1,m+1->s1,m+1
% for n = 2:N
%     for m = 2:M
%         rf(n,m)=(   rou001*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)) + ...
%                     rou011*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1) + ...
%                     rou101*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)) + ...
%                     rou111*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1)     ) ...
%                 ./...
%                 (   (1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)) + ...
%                     (1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1) + ...
%                     rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)) + ...
%                     rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1)     );
%     end
% end
% 
% % step6: s1,m->h1,m+1->s1,m+1
% for n = N:-1:2
%     for m = M:-1:2
%          rbl(n-1,m)=(   rou100*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         rou101*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)).*rbl(n,m).*rbr(n,m).*pi_in(n,m) + ...
%                         rou110*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         rou111*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1).*rbl(n,m).*rbr(n,m).*pi_in(n,m)     ) ...
%                     ./...
%                     (   (rou000+rou100)*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         (rou001+rou101)*(1- rf(n,m-1)).*(1- rbl(n,m-1)).*(1- pi_in(n,m-1)).*rbl(n,m).*rbr(n,m).*pi_in(n,m) + ...
%                         (rou010+rou110)*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         (rou011+rou111)*rf(n,m-1).*rbl(n,m-1).*pi_in(n,m-1).*rbl(n,m).*rbr(n,m).*pi_in(n,m)     );
% 
% 
%          rbr(n,m-1)=(   rou010*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         rou011*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*rbl(n,m).*rbr(n,m).*pi_in(n,m) + ...
%                         rou110*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         rou111*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*rbl(n,m).*rbr(n,m).*pi_in(n,m)     ) ...
%                     ./...
%                     (   (rou000+rou010)*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         (rou001+rou011)*(1- rf(n-1,m)).*(1- rbr(n-1,m)).*(1- pi_in(n-1,m)).*rbl(n,m).*rbr(n,m).*pi_in(n,m) + ...
%                         (rou100+rou110)*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*(1- rbl(n,m)).*(1- rbr(n,m)).*(1- pi_in(n,m)) + ...
%                         (rou101+rou111)*rf(n-1,m).*rbr(n-1,m).*pi_in(n-1,m).*rbl(n,m).*rbr(n,m).*pi_in(n,m)     );
%     end
% end
% if iter>iter_max
%     break;
% end
% iter = iter+1;
% end
% %% message passing: s->g->x->f
% pi_out = (rf.*rbr.*rbl)./((1.-rf).*(1.-rbr).*(1.-rbl)+(rf.*rbr.*rbl));
% 
% %%
% prob_S_temp = (pi_in.*pi_out)./((1.-pi_in).*(1.-pi_out)+(pi_in.*pi_out));
% prob_S = prob_S_temp(2:row+1, 2:col+1);
% 
% prob_S_reshape = reshape(prob_S,[row*col,1]);
% 
% % S_godec_1 = prob_S;
% % S_map = abs(reshape(S_godec_1,[row,col]));
% % figure(),colormap;imagesc(S_map);axis off;title ('prob_S')
% 
% ssd=1;
% 
% 
% 
% 
% 
% 
% 
% 
