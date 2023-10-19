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


H(1:I,:,:) = {zeros(R,T)};



%=====COM channels=====%
for i = 1:I
    for l = 1:C
        AOD2 = 2*rand()-1;
        AOA2 = 2*rand()-1;
        if l == 1
            AOA2 = AOD2;
        end
        AOD = AOD2*pi/2;
        AOA = AOA2*pi/2;
        AOAs = linspace(AOA-pi/18,AOA+pi/18, s_path);
        AODs = linspace(AOA-pi/18,AOA+pi/18, s_path);
        for s = 1:s_path
            strvctc1 = zeros(T,1);
            strvctc2 = zeros(d,1);
            for t = 1:T
                strvctc1(t) = exp(1i * (2*pi/0.0004) * (t-1) * 2 * sin(AOAs(s)));

            end
            for c = 1:d
                strvctc2(c) = exp(1i * (2*pi/0.0004) * (c-1) * 2 * sin(AODs(s)));
            end
            pow1 = trace(strvctc1*strvctc1');
            strvctc1 = strvctc1/sqrt(pow1);
            pow2 = trace(strvctc2*strvctc2');
            strvctc2 = strvctc2/sqrt(pow2);
            if l == 1 && s == 4
                xinta = sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
            end
            if l==1 && s~=4
                xinta = sqrt(0.1/2)*(randn(1,1)+1i*randn(1,1));
            end
            if l~=1 && s==4
                xinta = sqrt(0.3/2)*(randn(1,1)+1i*randn(1,1));
            end
            if l~=1 && s~=4
                xinta = sqrt(0.03/2)*(randn(1,1)+1i*randn(1,1));
            end

            Hc = xinta*strvctc2*strvctc1';
            H{i,K,K} = H{i,K,K}+Hc;
            
        end

    end

    H{i,K,K} = sqrt((T*R)/((s_path)*(C)))*H{i,K,K};

end


%=====clutter AOAS AODS=====%
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
%=====G_{\tau.s} and Gl=====%
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
            strvct1(t) = exp(1i * (2*pi/0.0003) * (t-1) * 2 * sin(LAODS(n)));
        end
        for s = 1:Nr
            strvct2(s) = exp(1i * (2*pi/0.0003) * (s-1) * 2 * sin(LAOAS(n)));
        end
        pow1 = trace(strvct1*strvct1');
        strvct1 = strvct1/sqrt(pow1);
        pow2 = trace(strvct2*strvct2');
        strvct2 = strvct2/sqrt(pow2);
        HR = strvct2*strvct1';
        if n==1
            sita = sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
            Hs1 = sita*HR;
            H{j,K,K} = Hs1;
        else
            sita = sqrt(10/2)*(randn(1,1)+1i*randn(1,1));
            Hl1 = sita*HR;
            Hl{j-I,n-1} = Hl1;
        end
    end


end




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