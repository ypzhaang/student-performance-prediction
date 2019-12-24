function [U,V,funs,tol_outs] = GRNMF_MM(W,M,U,V,dataGU,dataGM,para ,isNonNegative , type2Graph)
%% min_{U,V}  \| W.*(M  -U*V^T) \|_1 + \lambda_u \|U\|_F^2 + \lambda_v +\|V\|_F^2
%                      +alpha_u*Tr(U^T*L*U)+alpha_v*Tr(V^T*H*V)
%     s.t.    U>=0, V>=0          
% This script is for the papers:
% Zhouchen Lin, Chen Xu,and Hongbin Zha, Robust Matrix Factorization by Majorization-Minimization, T-PAMI, accepted
% Chen Xu, Zhouchen Lin,and Hongbin Zha, Relaxed Majorization-Minimization for Non-smooth and Non-convex Optimization, AAAI-16
% If you have any questions, feel free contract Chen Xu (xuen@pku.edu.cn)
%
% Copyright: Peking University
% Modifying by YUPEI ZHANG (@Northwestern Polytechnical University), Oct 3, 2018

%% Total Out, Inner Iterations and Tolerance
if isfield(para, 'max_out');      max_out = para.max_out;       else max_out = 300;   end
if isfield(para, 'max_in');       max_in  = para.max_in;        else max_in  = 5000;  end  
if isfield(para, 'tol_out');      tol_out = para.tol_out;       else tol_out = 1e-4;  end
if isfield(para, 'tol_in1');      tol_in1 = para.tol_in1;       else tol_in1 = 1e-5;  end
if isfield(para, 'tol_in2');      tol_in2 = para.tol_in2;       else tol_in2 = 1e-4;  end
if isfield(para, 'print_interval');      print_interval = para.print_interval;       else print_interval = 1;  end
 
%% Regularization Parameters
[m,n]=size(M);
if isfield(para, 'lambda_u');     lambda_u    = para.lambda_u;        else lambda_u   = 20/(m+n);   end
if isfield(para, 'lambda_v');     lambda_v    = para.lambda_v;        else lambda_v   = 20/(m+n);   end

%% Continuation by multiplying the upper bound parameter with Ratio until touch the upper bound 
if isfield(para, 'ratio');        ratio = para.ratio;                 else   ratio     = 0.004;     end
if isfield(para, 'scale_ratio');  scale_ratio = para.scale_ratio;     else scale_ratio = 1.2;       end
if isfield(para, 'max_ratio');    max_ratio = para.max_ratio;         else max_ratio   = 1;         end   
dLambda_u  = sum(W,2)';   
rdLambda_u = ratio*dLambda_u;
dLambda_v  = sum(W,1);  
rdLambda_v = ratio*dLambda_v;

%% Adaptive (rho) Penalty (Beta) for the inner LADMPSAP
if isfield(para, 'rho');        rho        = para.rho;                else rho       = 1.5;   end
if isfield(para, 'beta_type');  beta_type  = para.beta_type;          else beta_type = 2;     end
if isfield(para, 'beta_rate');  beta_rate  = para.beta_rate;          else beta_rate = 0.1;   end
Res = M-U*V';
nRes = norm(Res,'fro');
fun_res = sum(sum(abs(W.*Res)));
if beta_type == 0     % beta_rate =1e-3 is the defaut
   beta0 = beta_rate*fun_res/sum(W(:));
elseif beta_type == 1 % beta_rate =1e-1 is the defaut
   beta_start =  sum(abs(W.*Res),1)./sum(W,1);
   beta_start = sort(beta_start,'descend');
   beta_start = 1./beta_start;
   beta0 = beta_start(round(n* beta_rate))/3;    
elseif beta_type ==2  % beta_rate =1e-1 is the defaut
   beta0 = beta_rate*(m+n)*tol_in1;
end

%% Creating the graph L and H for U and V, respectively
opt_u=[];
opt_u.k=para.k_u;
opt_u.metric= type2Graph; %'HeatKernel';%'Cosine';'HeatKernel';
L=ConstructLap(dataGU, opt_u);
% L = full(L);

opt_v=[];
opt_v.k=para.k_u;
opt_v.metric= type2Graph; % 'HeatKernel';%'Cosine';
H=ConstructLap(dataGM, opt_v);
% H = full(H);

%% 需要调参的地方，非负的限制
funo =   + lambda_u*norm(U,'fro')^2 + lambda_v*norm(V,'fro')^2; %+ para.alpha_u * trace(U' * L * U) + para.alpha_v * trace(V' * H * V);
cW = 1 - W;
r = size(U,2);
iter_out = 1;
exceed_max_ratio =0;
funs = zeros(max_out,1);
tol_outs = zeros(max_out,1);
while 1
    if  iter_out == 1
        DeltaU = zeros(size(U));
        DeltaV = zeros(size(V));
        E = Res;
        Y = zeros(m,n);
    end 
    
    beta = beta0;
    U2norm = norm(U,2)^2;
    V2norm = norm(V,2)^2;
    Yhat = beta*(E+DeltaU*V'+U*DeltaV'-Res)+Y;
    for iter_in = 1:max_in 
        % Update E, DeltaU,and DeltaV In parallel
        temp = E-Yhat/(3*beta);  % 3 here is for 3 variables, refer to the papers for more details 
        E_temp = max(0,temp - 1/(3*beta)) + min(0,temp + 1/(3*beta));
        En = W.*E_temp+ cW.*temp;
        %% 这些需要调参
        DeltaUn = (-lambda_u*U+3*beta*V2norm*DeltaU-Yhat*V-para.alpha_u*(L*U))./((lambda_u+rdLambda_u+3*beta*V2norm)'*ones(1,r)+(para.alpha_u*L)*ones(m,r));
%         DeltaUn=max(0,DeltaUn);% Nonnegative
        DeltaVn = (-lambda_v*V+3*beta*U2norm*DeltaV-Yhat'*U-para.alpha_v*(H*V))./((lambda_v+rdLambda_v+3*beta*U2norm)'*ones(1,r)+(para.alpha_v*H)*ones(n,r));
%         DeltaVn=max(0,DeltaVn);% Nonnegative
        
        %% Relative Changes
        En_E = En-E; tnEn_E = En_E.*En_E; nEn_E = sqrt(sum(tnEn_E(:)));
        RelaChangeE = nEn_E/nRes;   %RelaChangeE = norm(E - En, 'fro')/ nRes;
        RelaChangeU = sqrt(V2norm)*norm(DeltaU-DeltaUn,'fro')/nRes;
        RelaChangeV = sqrt(U2norm)*norm(DeltaV-DeltaVn,'fro')/nRes;
        RelaChange = max([RelaChangeE,RelaChangeU,RelaChangeV]);
        
        E = En;
        DeltaU = DeltaUn;
        DeltaV = DeltaVn;
        
        %Update Y and Yhat
        Constaint_Error = E + DeltaU*V' + U*DeltaV'-Res;
        Y=Y + beta * Constaint_Error;
        Yhat = Y + beta * Constaint_Error ;
        if sqrt(beta) * RelaChange < tol_in1
            beta = min(1e5, rho * beta);
            tnCon = Constaint_Error .* Constaint_Error;   nCon = sqrt(sum(tnCon(:)));
            Rela_Constaint_Error = nCon/nRes;
            if Rela_Constaint_Error < tol_in2
                break
            end
        end
    end
    
   
   if isNonNegative == true
        U_n = max(0,U + DeltaU);
        V_n = max(0,V+ DeltaV);
   else
       if isNonNegative == false
           U_n = U + DeltaU;
           V_n = V + DeltaV;
       end
   end
  
    Res_n = M - U_n*V_n';
    fun = sum(sum(abs(W.*Res_n)))+lambda_u/2*norm(U_n,'fro')^2+lambda_v/2*norm(V_n,'fro')^2; %%+ para.alpha_u * trace(U_n' * L * U_n) + para.alpha_v * trace(V_n' * H * V_n);
    funs(iter_out,1) = fun;
    tol_outs(iter_out,1) = (funo-fun)/fun;
    %Ensuring the obj to monotonically decrease is enough for convergence from a trust region perspective
    if (funo-fun)/fun < tol_out  && (iter_out>= 10*(max(size(M))>800))  %The latter part is heuristic for the dinosaur tracks with nview =5   
        ratio=min(1,scale_ratio*ratio);                                
%         fprintf('Ratio increase to %1.2d \n', ratio)
        if ratio>= max_ratio
            exceed_max_ratio = exceed_max_ratio+1;
        end
        if exceed_max_ratio > 1
            disp('out_here');
            break
        end
        rdLambda_u = ratio*dLambda_u;
        rdLambda_v = ratio*dLambda_v;
    else
        U=U_n;
        V=V_n;
        Res= Res_n;
        nRes = norm(Res, 'fro');
        if mod(iter_out, print_interval)==0 || iter_out == 1
%            fprintf('iter_out=%d, fun = %1.3d, decrease = %1.2d \n',iter_out,fun, fun-funo)
        end
        funo=fun;
        if iter_out > max_out
           disp([' afer ',num2str(iter_out),' iters , converged']);
           break
        else
           iter_out = iter_out+1;
%            disp(['iter_ou now is : ', num2str(iter_out)]);
        end
    end
end