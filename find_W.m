function W = find_W(U,H,V,T,Nr,I,S,K,d)
W = cell(I+S,K);

W(1:I+S) = {zeros()};

for i=1:I
    for k=1:K
        W{i,k} = inv(eye(d)-U{i,k}'*H{i,k,k}*V{k,i});

    end
end
for i=I+1:I+S
    for k=1:K
    W{i,k} = inv(eye(Nr/S)-U{i,k}'*H{i,k,k}*V{k,i});

    end
end
end