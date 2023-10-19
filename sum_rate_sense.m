function [sense_rate] = sum_rate_sense(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,alpha1)
rate = zeros(I+S,K);
for i=I+1:I+S
    for k = 1:K
        tempss = zeros();
        for l=1:I+S
            for j=1:K
                if l ~= i

                    tempss = tempss + H{i,K,j}*(V{j,l}*V{j,l}')*(H{i,K,j}');
                end
            end
        end
        % in R_{\tau}, com signal in G_{\tau,s} was ignored
%         for m = I+1 : I+S
%             if m ~= i
%                 for t = 1:I+S
%                 tempss = tempss + H{m,k,j}*(V{K,t}*V{K,t}')*H{m,k,j}';
%                 end
%             end
%         end
    end
    tempct = zeros();

        for t = 1 : L
            for u = 1:I+S
                tempct = tempct + Hl{1,t}*(V{K,u}*V{K,u}')*Hl{1,t}';
            end
        end



    rate(i,K) = log2(det(eye(Nr)+H{i,K,K}*(V{K,i}*V{K,i}')*H{i,K,K}'*( ...
        sigma2*eye(Nr)+tempct)^-1));

end


sense_rate = real(sum(rate(I+1:I+S)))/S;

end