function [rate_saver_xls,average_rate_xls,sense_rate_xls,com_rate_xls,MMSE_rate_xls, MMSE_average_xls,MMSE_sense_xls,MMSE_com_xls] = init(snr)

K = 1; % 基站个数
T = 16; % 发射天线个数
R = 4; % 接收天线个数
d = R;
Nr = 4;
epsilon = 1e-3; % 收敛条件
sigma2 =1; % 噪声功率
% snr = 25; % 信噪比
P = db2pow(snr)*sigma2; % 发射功率
SINRZF = 25;
I = 3; % 每个基站服务的用户个数
S = 1;%感知目标个数
L = 3;%感知目标的杂波(clutter)个数
s_path = 9;
C = 1;
alpha1 = ones(I+S,K); % 权重系数，都假设相同
% a = sqrt((T*Nr)/(L+1));
count = 0;

 % 假设每个用户有R路独立的数据流

max_iter = 100;
rate = []; % 初始化一个空向量记录rate

H = cell(I+S,K,K); % 信道系数
Hl = cell(S,L);
Hs = cell(I+S,K,K);

H(1:I,:,:) = {zeros(R,T)};




for i = 1:I
    for l = 1:C 
        for s = 1:s_path+1
        AOD2 = 2*rand()-1;
        AOA2 = 2*rand()-1;
        if l == 1
            AOA2 = AOD2;
        end
        AOD = AOD2*pi/2;
        AOA = AOA2*pi/2;
       
       
            strvctc1 = zeros(T,1);
            strvctc2 = zeros(d,1);
            for t = 1:T
                strvctc1(t) = exp( -1i* (2*pi*0.5) * (t-1) * sin(AOA));

            end
            for c = 1:d
                strvctc2(c) = exp(-1i * (2*pi*0.5) * (c-1) * sin(AOD));
            end
            pow1 = trace(strvctc1*strvctc1');
            strvctc1 = strvctc1/sqrt(pow1);
            pow2 = trace(strvctc2*strvctc2');
            strvctc2 = strvctc2/sqrt(pow2);
            if l == 1 && s == 4
                xinta = sqrt(1/2)*(normrnd(0,1)+1i*normrnd(0,1));
            end
            if l==1 && s~=4
                xinta = sqrt(1/2)*(normrnd(0,0.1)+1i*normrnd(0,0.1));
            end
            if l~=1 && s==4
                xinta = sqrt(1/2)*(normrnd(0,0.1)+1i*normrnd(0,0.1));
            end
            if l~=1 && s~=4
                xinta = sqrt(1/2)*(normrnd(0,0.1)+1i*normrnd(0,0.1));
            end

            Hc = xinta*strvctc2*strvctc1';
            H{i,K,K} = H{i,K,K}+Hc;
            
        end

    end

    H{i,K,K} = sqrt((T*R)/((s_path+1)*(C)))*H{i,K,K};

end


Hs(1:I+S,:,:) = {zeros(Nr,T)};
Hl(:,:) = {zeros(Nr,T)};
LAOA = [];
LAOD = [];
for cl = 1:L
    AOAL = 2*rand()-1;
    AODL = 2*rand()-1;
    AODL = AODL*pi/2;
    AOAL = AOAL*pi/2;
    LAOA = [LAOA AOAL];
    LAOD = [LAOD AODL];
end
for j = I+1 : I+S

       
       
       LAOAS = LAOA;
       LAODS = LAOD;
        % % % % % % % % % % %
        AOD0 = 2*rand()-1;
        AOA0 = 2*rand()-1;

        AOD = AOD0*pi/2;
        AOA = AOA0*pi/2;
        LAOAS = [AOA LAOA];
        LAODS = [AOD LAOD];
     for n = 1:L+1
        strvct1 = zeros(T,1);
        strvct2 = zeros(Nr,1);
        for t = 1:T
            strvct1(t) = exp(-1i * (2*pi*0.5) * (t-1) * sin(LAODS(n)));
        end
        for s = 1:Nr
            strvct2(s) = exp(-1i * (2*pi*0.5) * (s-1) * sin(LAOAS(n)));
        end
        pow1 = trace(strvct1*strvct1');
        strvct1 = strvct1/sqrt(pow1);
        pow2 = trace(strvct2*strvct2');
        strvct2 = strvct2/sqrt(pow2);
        HR = strvct2*strvct1';
        if n==1
            sita = sqrt(1/2)*(normrnd(0,1)+1i*normrnd(0,1));
            Hs1 = sita*HR;
            H{j,K,K} = Hs1;
        else
            sita = sqrt(1/2)*(normrnd(0,1)+1i*normrnd(0,1));
            Hl1 = sita*HR;
            Hl{j-I,n-1} = Hl1;
        end
    end


end




U = cell(I+S,K);
U(1:I)={zeros()};
U(I+1,:) = {zeros()};
% ===============================================================================
% 随机初始化发射波束向量

V_ini = V_init_v2(I,S,K,T,Nr,d,P);
V = V_ini;



% 求初始化发射波束V后求系统和速率
rate_old = sum_rate_all(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
average_rate_old = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
rate = [];
average_rate = [];
sense_rate = [];
% sense_rate_old = sum_rate_sense(H,V,Hzf,Vzf,sigma2,T,R,I,K,alpha1);
% sense_rate = [sense_rate sense_rate_old];
com_rate = [];
% com_rate_old = sum_rate_com(H,V,Hzf,Vzf,sigma2,T,R,I,K,alpha1);
% com_rate = [com_rate com_rate_old];


sense_rate_saver = [];
sense_rate_saver_zf = [];
com_rate_saver = [];
trE_saver = [];
com_rate_saver_zf = [];
rate_saver = [];
average_rate_saver = [];
x_axis = [];

% ========================================================================

% for w = 0.00:0.01:0.99

    w = 0.99;
    wr = w;
    wc = (1-wr);

    %     Vzf = mua * HzfM * ((HzfM'*HzfM)^-1);
    %     Vzf(:,1:I*R) = wc/I*Vzf(:,1:I*R);
    %     Vzf(:,I*R+1) = wr*Vzf(:,I*R+1);
    %     Vzf = mat2cell(Vzf,T,[R R 1]);


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
    %     sense_rate_zf = [];
    sense_rate_old = sum_rate_sense(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    sense_rate = [sense_rate sense_rate_old];
    %     sense_rate_zf = [sense_rate_zf sense_rate_zf_old];
    com_rate = [];
    %     com_rate_zf = [];
    com_rate_old = sum_rate_com(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
    com_rate = [com_rate com_rate_old];
    %     com_rate_zf = [com_rate_zf,com_rate_zf_old];
    %     alpha1(I+1,K) = wr;
    iter1 = 1;
    while(1)

        U = find_U(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,d); % Tbale I line 4 in p.4435
        W = find_W(U,H,V,T,Nr,I,S,K,d); % Tbale I line 5 in p.4435

        V = find_V(alpha1,H,Hl,U,W,T,R,I,S,L,K,P); % Tbale I line 6 in p.4435
        rate_new = sum_rate_all(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        average_rate_new = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        sense_rate_new = sum_rate_sense(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        sense_rate = [sense_rate sense_rate_new];
        %         sense_rate_zf = [sense_rate_zf sense_rate_zf_new];
        com_rate_new= sum_rate_com(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1);
        com_rate = [com_rate com_rate_new];
        %         com_rate_zf = [com_rate_zf com_rate_zf_new];

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
   
rate_saver_xls = reshape(rate_saver,[1,1]);
sense_rate_xls = reshape(sense_rate_saver,[1,1]);
com_rate_xls = reshape(com_rate_saver,[1,1]);
average_rate_xls = reshape(average_rate_saver,[1,1]);
% 
[MMSE_rate,MMSE_average,MMSE_com,MMSE_sense] = MMSE_baseline(H,Hl,snr,V_ini);
MMSE_rate_xls = reshape(MMSE_rate,[1,1]);
MMSE_average_xls = reshape(MMSE_average,[1,1]);
MMSE_com_xls = reshape(MMSE_com,[1,1]);
MMSE_sense_xls = reshape(MMSE_sense,[1,1]);
