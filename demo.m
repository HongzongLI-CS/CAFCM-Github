addpath kit
%% load dataset
id=1; %id=1, id=2, ..., id=14 denote different datasets
[X,label,p] = load_dataset(id);
n=size(X,1);
m=size(X,2);
% n is the number of data points, 
% m is the number of features,
% p is the number of clusters. 

%% hyper-parameters setting
M=5;
N=2;

%% cooling schedule selection
schedule=2; 
% 1:exponential cooling schedule 
% 2:polynomial cooling schedule

[resulting_label,obj_value,time_cost]=cafcm(X,n,m,p,schedule,M,N);
%% It is also optional to use default parameters
%[resulting_label,obj_value,time_cost]=cafcm(X,n,m,p,schedule);
%[resulting_label,obj_value,time_cost]=cafcm(X,n,m,p);

%% Save results
%pre = [cd,'/results/'];
filename = ['/example_',num2str(id),'_cafcm.txt'];
savePath = [cd,filename];
writematrix([obj_value,time_cost,resulting_label'],savePath,'Delimiter','\t','WriteMode','append')