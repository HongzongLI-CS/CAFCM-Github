function [X,label,p] = load_dataset(id)
addpath datasets
%% load data
if id==1
    data=load('winequality-red.dat');
    p=6;
    X=data(:,1:end-1);
    label=data(:,end)';
    %optimum_value=
elseif id==2
    data=load('segment.dat');
    p=7;
    X=data(:,1:end-1);
    X(:,3)=[];
    label=data(:,end)';
elseif id==3
    data=load('spambase.dat');
    p=2;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==4
    data=load('winequality-white.dat');
    p=11;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==5
    data=load('banana.dat');
    p=2;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==6
    data=load('phoneme.dat');
    p=2;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==7
    data=load('page-blocks.dat');
    p=5;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==8
    data=load('satimage.dat');
    p=7;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==9
    data=load('coil2000.dat');
    p=2;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==10
    data=load('penbased.dat');
    p=10;
    X=data(:,1:end-1);
    label=data(:,end)';
elseif id==11
    data=load('EGS.csv');
    p=2;
    X=data(:,1:end-1);
    label=data(:,end)';
    label=label+1;
elseif id==12
    data=load('warpPIE10P.mat');
    p=10;
    X=data.X;
    label=data.Y';
elseif id==13
    data=load('nci9.mat');
    p=9;
    X=data.X;
    label=data.Y';
elseif id==14
    data=load('optdigits.dat');
    p=10;
    X=data(:,1:end-1);
    X(:,1)=[];
    X(:,39)=[];
    label=data(:,end)';
end

X=(X-min(X))./(max(X)-min(X));

end

