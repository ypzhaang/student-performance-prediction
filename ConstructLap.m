function Lap = ConstructLap(Y,opt)
%CONSTRUCTLAP Summary of this function goes here
%   Detailed explanation goes here
% Y: the data with one row being a point
if isempty(Y)  
    warning('Y matrix is empty!'); 
    return; 
end

if isfield(opt,'K') 
    opt.k=5; 
end

options=[];
options.k=opt.k;
options.WeightMode=opt.metric;
W=constructW(Y,options);
W=(W+W')/2;

% Lap  = W;
Lap=spdiags(full(sum(W,2)),0,size(Y,1),size(Y,1))-W;

end

