clear;
[trainMat , testMat , dataTest , dataTrain , row , col ,dataGU , dataGM] = MovieDataSmallTest();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
NMFErr = [];
MFrmse = [];
MFErr = [];
GRMFrmse = [];
GRMFErr = [];
     %% GRMF -- when the                                                
for alpha_v = 0.24 : 0.02 :1.2 % the values of all the parameters could changed as you like
      [ M_mm ,Rmse,Acc ,Err1 , funs, tol_outs] = runGRMF(dataTest,trainMat ,testMat , dataGU, dataGM , 5 , 15 ,2.74, 0.50 , 0.1, 300 , 1000 , false ,'Cosine');
      disp([' when use the GRMF, alpha_u =  ' ,num2str(0.5) ,' alpha_v =  ',num2str(0.54),' RMSE =  ' num2str(Rmse) , ' acc =  ' , num2str(Acc) , ' Err1 =  ' , num2str(Err1)]);
      GRMFrmse = [GRMFrmse ; Rmse];
      GRMFErr = [GRMFErr ; Err1]; 
      plot(funs);plot on; plot(tol_outs);
      xlabel('iter times');
      ylabel('Objective Value');
      
      % RMF
      [ M_mm ,Rmse,Acc ,Err1 , funs, tol_outs] = runGRMF(dataTest,trainMat ,testMat , dataGU, dataGM , 4 , 0 ,3.9, 0 ,0  ,30 , 1250 , false ,'HeatKernel');
      disp([' when use the  RMF ', 'RMSE = ' num2str(Rmse) , ' acc = ' , num2str(Acc) ,  ' Err1 = ' , num2str(Err1)]);
      RMFrmse = [RMFrmse ; Rmse];
      RMFErr = [RMFErr ; Err1];

      % GNMF                              
      [ GNMF_all , Rmse , Acc , Err1] = runGNMF( dataTest ,trainMat , testMat , row , col , 3000 , 100 , 5);
      disp([' when use the GNMF the paras is the best,then  GNMF ' ,' RMSE = ' num2str(Rmse) , ' acc = ' , num2str(Acc) ,  ' Err1 = ' , num2str(Err1)]);
      GNMFrmse = [GNMFrmse ; Rmse];
      GNMFErr = [GNMFErr ; Err1];

      % MF
      [ MF_all , Rmse , Acc , Err1] = runMF( dataTest ,trainMat , testMat , row , col , 3000  , 5);
      disp([' when use the MF nClass =  ' ,' RMSE = ' num2str(Rmse) , ' acc = ' , num2str(Acc) ,  ' Err1 = ' , num2str(Err1)]);
      MFrmse = [MFrmse ; Rmse];
      MFErr = [MFErr ; Err1];
end

%% Visualization after running the algorithm
    % GRALSErr.mat & nclass-GRMFErr.mat could be collected by yourself
    hold on;
    x = 3:13;

    GRALSErr = open('GRALSErr.mat');
    GRALSErr = GRALSErr.GRALSErr;
    GRMFErr = open('nclass-GRMFErr.mat');
    GRMFErr = GRMFErr.GRMFErr;
    plot(x , GRMFErr(1:11)  , 'r-o' , 'LineWidth' , 2);
    plot(x , GRALSErr(1:11) , 'b-*' , 'LineWidth' , 2);
    plot(x , RMFErr(1:11)   , 'g-d' , 'LineWidth' , 2);
%     plot(x , GNMFErr , 'c-x');
%     plot(x , RMFErr(1:11)  , 'k-s');

    xlabel('rank');
    ylabel('Err1');

    legend( 'GRMF' ,'GRALS' , 'RMF'  , 'Location' , 'northeast');
    
%% Visualization of the Fig 3.(a)
    mesh(GRMFErr)
    colormap copper
    colormap bone
    colormap cool
    colormap hsv
    zlabel('Err1')
    xlabel('alpha_u')
    xlabel('alpha-u')
    ylabel('alpha-v')