function [outMat] = MF(trainMat, k ,steps,alpha,beta )
%MATRIXFACTORIZATION ´Ë
[row,col] = size(trainMat);
P = rand(row , k) * 10;
Q = rand(col , k) * 10;
Q = Q';
for step = 1 : steps
    for i = 1 : row
        for j = 1 : col
            if trainMat(i,j) > 0
                eij = trainMat(i,j) - P(i , :) * Q(:,j);
                for kit = 1 : k
                    P(i,kit) = P(i,kit) + alpha * (2 * eij * Q(kit,j) - beta * P(i,kit));
                    Q(kit,j) = Q(kit,j) + alpha * (2 * eij * P(i,kit) - beta * Q(kit,j));
                end
            end
        end
    end
    outMat = P * Q;
    e = 0;
    for i = 1 : row
        for j = 1 : col
            if trainMat(i,j) > 0 
                e = e + (trainMat(i,j) - P(i,:) * Q(:,j))^2;
                for kit = 1 : k
                    e = e + (beta / 2) * (P(i,kit)^2 + Q(kit,j)^2);
                end
            end
        end
    end
    if e < 0.001
        break
    end
end
end

