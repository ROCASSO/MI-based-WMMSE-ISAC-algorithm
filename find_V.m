function V = find_V(alpha1,H,Hl,U,W,T,R,I,S,L,K,P)
    V = cell(K,I+S);
    A = cell(K);
    B = cell(K);
    
    A(:) = {zeros()};
    B(:) = {zeros()};
    
    
    for k=1:K     
        for j=1:K
            for l=1:I+S
                
                A{k} = A{k} + alpha1(l,j)*H{l,j,k}'*U{l,j}*W{l,j}*(U{l,j}')*H{l,j,k};
            end
        end   
    end 
   
    for c = I+1 : I+S
        s = c-I;
        for l = 1:L
            B{K} = B{K} + alpha1(c,K)*(Hl{s,l}')*U{c,K}*W{c,K}*(U{c,K}')*Hl{s,l};
        end
    end



        
    
    max_iter = 1000; % Bisection
    mu = zeros(K,1);
    for k=1:K 
        mu_min = 0;
        mu_max = 1e7;
        iter = 0;
        while(1)
            mu1 = (mu_max+mu_min) / 2;
            P_tem = 0;
            for i=1:I+S 
                V_tem = inv((A{k}+B{k}+mu1*eye(T)))*(alpha1(i,k)*((H{i,k,k}')*U{i,k}*W{i,k}) ...
                    );
                P_tem = P_tem + real(trace(V_tem*V_tem'));
            end
            if P_tem > P
                mu_min = mu1;
            else
                mu_max = mu1;
            end
            iter = iter + 1;

            if abs(mu_max - mu_min) < 1e-7 || iter > max_iter
                break
            end
        end
        mu(k) = mu1;
    end

    for i=1:I+S

        for k=1:K
            V{K,i} =  inv((A{k}+B{k}+mu1*eye(T)))*(alpha1(i,k)*((H{i,k,k}')*U{i,k}*W{i,k}) ...
                      );
                 

        end 
    end

end