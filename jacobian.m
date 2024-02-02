function J = jacobian(a, L,alpha)
m = size(a,1) -1;
J = zeros(m+1,m+1);
for j = 0:m
    for k= 0:m
        mu =0;
        sum_else = 0;
        for k2 = -m:m
            for k3 = -m:m
                index_2 = abs(k2);
                index_3 = abs(k3);
                if k2+k3 == k-j || k2+k3== k+j
                    sum_else = sum_else +a(index_2+1)*a(index_3+1);

                end

            end
        end
        if k == j 
            mu = -(((k)*pi)/L)^4 - 2*((k*pi)/L)^2 - (1-alpha);

        end 
        J(k+1,j+1) =mu -3*sum_else;
        

    end
end

return