function [rate_saver_xls,average_rate_xls,sense_rate_xls,com_rate_xls] = algorithm(snr)


 

K = 1; %Number of BSt (for future works if possible)
T = 16; % Nt
R = 4; % UE antenna
d = R;
Nr = 4;%APr antenna
epsilon = 1e-3; % tolerance
sigma2 =1; % noise
snr = 25; % 信噪比
P = db2pow(snr)*sigma2; % power
SINRZF = 25;
I = 3; % UE
S = 1;%ST
L = 3;%clutter patches
s_path = 10;%paths in each clusters
C = 4;%clusters
alpha1 = ones(I+S,K); % weights





max_iter = 30;
rate = []; 

H = cell(I+S,K,K); 
Hl = cell(S,L);%clutter channel(response matrix)
%set the channel based on SV-model
U = cell(I+S,K);
U(1:I)={zeros()};
U(I+1,:) = {zeros()};
% ===============================================================================


V_ini = V_init_v2(I,S,K,T,Nr,d,P);
V = V_ini;




rate_old = sum_rate_all(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
average_rate_old = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
rate = [];
average_rate = [];
sense_rate = [];

com_rate = [];



sense_rate_saver = [];
sense_rate_saver_zf = [];
com_rate_saver = [];
com_rate_saver_zf = [];
rate_saver = [];
average_rate_saver = [];
x_axis = [];
% ========================================================================

for w = 0.00:0.01:0.99


    wr = w
    wc = (1-wr);


% for the following I+1:I+S, it was prepared for further investigation
% about multiple STs if possible
    V = V_ini;
    alpha1(I+1:I+S,K) = wr/S;
    for j=1:K
        for l=1:I
            alpha1(l,j) = wc/(I);
        end
    end
    rate_old = sum_rate_all(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    average_rate_old = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    rate = [rate rate_old];
    average_rate = [average_rate average_rate_old];
    sense_rate = [];
    
    sense_rate_old = sum_rate_sense(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    sense_rate = [sense_rate sense_rate_old];
    
    com_rate = [];
    
    com_rate_old = sum_rate_com(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    com_rate = [com_rate com_rate_old];
    
    iter1 = 1;
    while(1)
       %algorithm
        U = find_U(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,d); 
        W = find_W(U,H,V,T,Nr,I,S,K,d); 
        V = find_V(alpha1,H,Hl,U,W,T,R,I,S,L,K,P);
        rate_new = sum_rate_all(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        average_rate_new = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        sense_rate_new = sum_rate_sense(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        sense_rate = [sense_rate sense_rate_new];
       
        com_rate_new= sum_rate_com(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        com_rate = [com_rate com_rate_new];
    

        iter1 = iter1 + 1;
        if abs(rate_new-rate_old) / rate_old < epsilon|| iter1 > max_iter

            break;
        end
        rate_old = rate_new;
        sense_rate_old = sense_rate_new;
        com_rate_old = com_rate_new;
        average_rate_old = average_rate_new;
    end
    rate_saver = [rate_saver rate_new];
    average_rate_saver = [average_rate_saver average_rate_new];

    sense_rate_saver = [sense_rate_saver sense_rate_new];


    com_rate_saver = [com_rate_saver com_rate_new];

    
    x_axis = [x_axis wr];
    
end
rate_saver_xls = reshape(rate_saver,[100,1]);
sense_rate_xls = reshape(sense_rate_saver,[100,1]);
com_rate_xls = reshape(com_rate_saver,[100,1]);
average_rate_xls = reshape(average_rate_saver,[100,1]);
