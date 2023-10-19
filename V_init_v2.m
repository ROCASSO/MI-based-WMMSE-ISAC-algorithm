function [V_ini] = V_init_v2(I,S,K,T,Nr,d,P)
V = cell(K,I+S); 


for i=1:I
    for k=1:K
        v = normrnd(0,1,T,d)+1i*normrnd(0,1,T,d); 
        V{K,i}=v;
    end
end
for i=I+1:I+S
    for k=1:K
        v = normrnd(0,1,T,Nr)+1i*normrnd(0,1,T,Nr); 
        V{K,i}=v;
    end
end

B = cell2mat(V);
pw = trace(B*B');
for i = 1:I+S

pw = trace(V{K,i}*V{K,i}');
V{K,i} = sqrt(P/I+S)*sqrt(pw)^-1*V{K,i};
end
V_ini = V;
end