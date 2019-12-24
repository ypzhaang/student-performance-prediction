function [ M_mm , Rmse , Acc , Err1] = runGNMF( dataTest ,trainMat , testMat, row , col , maxIter , alpha ,nClass)
%RUN 此处显示有关此函数的摘要
%   此处显示详细说明
%GNMF learning
fea = trainMat;
% fea = NormalizeFea(fea);
options = [];
options.WeightMode = 'HeatKernel';  
W = constructW(fea,options);
options.maxIter = maxIter;%100;
options.alpha = alpha ;% 100;
rand('twister',5489);
[U,V] = GNMF(fea',nClass,W,options); %'

%Clustering in the GNMF subspace
M_mm = V * U';

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

