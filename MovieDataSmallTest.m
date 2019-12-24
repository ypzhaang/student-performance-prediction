function [trainMat , testMat , dataTest , dataTrain , row , col ,dataGU , dataGM] = MovieDataSmallTest()
%% 
clear;

%% load data
data = load('u.data');
data = data(:,1:3);
data = data(data(: , 1) < 200 & data(: , 2) < 300,:);

dataGU = load('userMat.data');
dataGU = dataGU';
[rowG , colG] = size(dataGU);
dataGU = dataGU(dataGU(:,1) < 200 , :);
dataGU = dataGU(: , 2 : colG);


dataGM = load('movieInfo.data'); 
[rowM , colM] = size(dataGM);
dataGM = dataGM(dataGM(:,1) < 300 , :);

dataGM = dataGM(:,2 : colM);


%% get train and test
len = length(data);

index = randperm(len);
rand('twister',5489);
shit = data;

for i = 1 : len
    data(i,:) = shit(index(i),:);
end

ratio = 0.1;

dataTest = data(1 : int32(len * ratio) , :);
dataTrain = data(int32(len * ratio) + 1 : len , :);

%% initilization testMat and TrainMat

row = length(unique(data(:,1)));
col = length(unique(data(:,2)));

trainMat = zeros(row , col);
testMat = zeros(row , col);

for i = 1 : length(dataTrain)
    uid = dataTrain(i,1);
    mid = dataTrain(i,2);
    rate = dataTrain(i,3);
    trainMat(uid , mid) = rate;
end

for i = 1 : length(dataTest)
    uid = dataTest(i,1);
    mid = dataTest(i,2);
    rate = dataTest(i,3);
    testMat(uid , mid) = rate;
end
end

