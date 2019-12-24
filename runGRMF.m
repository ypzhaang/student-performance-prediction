function [ M_mm ,Rmse,Acc ,Err1 , funs, tol_outs] = runGRMF(dataTest,trainMat ,testMat , dataGU, dataGM , nClass , k ,lambda, alpha_u , alpha_v  ,max_out , max_in , isNonNegative , type2Graph)
%     nClass =11;
    M = trainMat;
    [row,col] = size(trainMat);
    W = ones(row,col);
    ze_ro = find(M == 0);
    W(ze_ro) = 0;
 
    M_cen = M;
    M_cen = W.*M_cen;

    % Inilialization
    [U,S,V] = mySVD(M_cen,nClass); 
    U0 = U*sqrt(S);
    V0 = V*sqrt(S);

    clear para
    [m,n] = size(M);
%     lambda = 1020;%620;
    para.lambda_u = (lambda);
    para.lambda_v = (lambda);
%     alpha = 0.55;
    para.alpha_u = alpha_u;
    para.alpha_v = alpha_v;
%     k = 9; 
    para.k_u = k;
    para.k_v = k;
    para.rho = 3.2;
    para.max_out = max_out;
    para.max_in = max_in;
    % disp(tic) %% 
    tic
    [U_mm, V_mm, funs, tol_outs] = GRMF_MM(W , M_cen , U0 , V0 , dataGU , dataGM , para ,isNonNegative , type2Graph); %% 加入DataG 作出一个新图
    
    if(length(find(U_mm < 0) ) > 1 && isNonNegative == true)
        disp('there are wrongs with U_mm Negative.');
    end
    if(length(find(V_mm < 0) ) > 1 && isNonNegative == true)
        disp('there are wrongs with V_mm Negative.');
    end
    %% 评估
    M_mm = U_mm * V_mm';
    
    if length(dataTest) ~= 1
        outMat = zeros(row , col);
        W = zeros(row , col);
        for i = 1 : length(dataTest)
            uid = dataTest(i,1);
            mid = dataTest(i,2);
            outMat(uid , mid) = M_mm(uid , mid);
            W(uid , mid ) = 1;
        end
    end
    if length(testMat) ~= 1
        [Rmse , Acc] = OrderRmseAcc(outMat , testMat);
        Err1 = sum(sum(abs(W.*(outMat-testMat))))/sum(W(:));
    end
end

