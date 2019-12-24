function [ M_mm , Rmse , Acc , Err1] = runMF( dataTest ,trainMat , testMat, row , col , maxIter  ,nClass)
%RUNMF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

k = nClass;
steps = maxIter;
alpha = 0.001;
beta = 0.02;
M_mm = MF(trainMat, k ,steps,alpha,beta );

% [M_mm] = MF(trainMat, nClass ,maxIter,alpha,beta );
outMat = zeros(row , col);
W = zeros(row , col);
for i = 1 : length(dataTest)
    uid = dataTest(i,1);
    mid = dataTest(i,2);
    outMat(uid , mid) = M_mm(uid , mid);
    W(uid , mid ) = 1;
end
[Rmse , Acc] = OrderRmseAcc(outMat , testMat);
Err1 = sum(sum(abs(W.*(outMat-testMat))))/sum(W(:));
end

