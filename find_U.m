function U = find_U(H,Hl,V,sigma2,T,Nr,R,I,S,L,K,d)
J = cell(I+S,K);
J(:)={zeros()};
U = cell(I+S,K);
U(1:I+S)={zeros()};

for i=1:I
    for k=1:K
        for j=1:K
            for l=1:I+S
                J{i,k} = J{i,k} + H{i,K,K}*V{K,l}*(V{K,l}')*(H{i,K,j}'); % 算法Table I, 第四行括号求和的部分
            end
        end
        J{i,k} = J{i,k} + sigma2*eye(R); 
        U{i,k} = J{i,k}\H{i,k,k}*V{k,i}; 
    end
end

for c = I+1:I+S
    s = c-I;
    for k=1:K
        tempss = zeros();

        for j = 1:I+S

            tempss = tempss + H{c,K,K}*V{K,j}*(V{K,j}')*(H{c,K,K}');
        end

        temp = zeros();
        for t = 1:L

            for l=1:I+S
                temp = temp + Hl{s,t}*V{K,l}*V{K,l}'*Hl{s,t}';

            end
        end
% 
        J{c,k} = temp+tempss;
        J{c,k} = J{c,K} + sigma2*eye(Nr); 
        U{c,k} = J{c,K}\H{c,K,K}*V{K,c}; 
    end
end



end