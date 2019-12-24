function [rmse,acc] = OrderRmseAcc(outMat,testMat)
%ORDERRMSEACC 此处显示有关此函数的摘要
%   此处显示详细说明
[rowT,colT] = size(testMat);
[rowO,colO] = size(outMat);
if [rowT,colT] == [rowO,colO]
    sumRmse = 0;
    numRmse = 0;
    numAcc = 0;
    %可以继续
    for i = 1 : rowT
        for j = 1 : colT
            if testMat(i,j) ~= 0
                mid = (testMat(i,j) - outMat(i,j))^2;
                sumRmse = sumRmse + mid;
                numRmse = numRmse + 1;  
                
                if int32(testMat(i,j)) == int32(outMat(i,j))
%                 if abs(testMat(i,j) - outMat(i,j)) <= 0.2
                    numAcc = numAcc + 1;
                end
            end
        end
    end
    acc = numAcc / numRmse;
    rmse = sqrt(sumRmse / numRmse);
%     disp(num2str(numRmse));
else
    disp(['testMat 的大小为 = ', [rowT,colT]]);
    disp(['outMat 的大小为 = ', [rowO,colO]]);
end

end

