function [rmse,acc] = OrderRmseAcc(outMat,testMat)
%ORDERRMSEACC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[rowT,colT] = size(testMat);
[rowO,colO] = size(outMat);
if [rowT,colT] == [rowO,colO]
    sumRmse = 0;
    numRmse = 0;
    numAcc = 0;
    %���Լ���
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
    disp(['testMat �Ĵ�СΪ = ', [rowT,colT]]);
    disp(['outMat �Ĵ�СΪ = ', [rowO,colO]]);
end

end

