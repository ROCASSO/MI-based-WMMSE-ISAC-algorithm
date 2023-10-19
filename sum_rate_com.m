function [com_rate] = sum_rate_com(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1)
%average CMI per user%
rate = zeros(I+S,K);

for i=1:I
    for k = 1:K
        temp = zeros();
        for l=1:I+S
            for j=1:K
                if l~=i || j~=k
                    temp = temp + H{i,k,j}*V{j,l}*(V{j,l}')*(H{i,k,j}');
                end
            end
        end
        rate(i,k) = log2(det(eye(R)+H{i,k,k}*V{k,i}*(V{k,i}')*(H{i,k,k}') ...
            *inv(temp + sigma2*eye(R))));
        
    end
end
com_rate = real(sum(rate(1:I)))/I;
end