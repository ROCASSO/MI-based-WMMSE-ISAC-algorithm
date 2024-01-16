%WMMSE-ISAC Mutual Information-Based Integrated Sensing and
% Communications: A WMMSE Framework
% Yizhou Peng, Songjie Yang, Wanting Lyu, Ya Li, Hongjun He,
% Zhongpei Zhang, Member, IEEE, and Chadi Assi, Fellow, IEEE
%randomly chosen your seed
sample_num = 300;
x_axis_all = 0:0.01:0.99;
rate_all = zeros(100,1,sample_num);
com_rate_all = zeros(100,1,sample_num);
sense_rate_all = zeros(100,1,sample_num);
samples = zeros(1,100,sample_num);
for snr = 5:5:35
    for sampletime = 1:sample_num
        [R,S,C] = algorithm(snr);

        rate_all(:,:,sampletime) = R;
        com_rate_all(:,:,sampletime) = C;
        sense_rate_all(:,:,sampletime) = S;
       

    end



rate_average = sum(rate_all,3)/sample_num;
com_rate_average = sum(com_rate_all,3)/sample_num;
sense_rate_average = sum(sense_rate_all,3)/sample_num;
end
