function [average_system_rate] = sum_rate(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1)
%=====average CMI per user + SMI=====%
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
        rate(i,k) = log2(det(eye(R)+(H{i,k,k}*V{k,i}*(V{k,i}')*(H{i,k,k}') ...
            *inv(temp + sigma2*eye(R)))));
        
    end
end
for i=I+1:I+S
  

        tempct = zeros();
        for t = 1 : L
            for u = 1:I+S
                tempct = tempct + Hl{1,t}*V{K,u}*V{K,u}'*Hl{1,t}';
            end
        end



    rate(i,K) = log2(det(eye(Nr)+H{i,K,K}*V{K,i}*(V{K,i})'*(H{i,K,K})'*inv( ...
        sigma2*eye(Nr) +tempct)));
end
sense_rate = real(sum(rate(I+1:I+S)))/S;
com_rate = real(sum(rate(1:I)))/I;
average_system_rate = com_rate+sense_rate;
end