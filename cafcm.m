function [label,gbest,all_time]=cafcm(X,n,m,p,schedule,M,N)
if nargin<7
    % default parameters
    M=5;
    N=2;
end
if nargin<5
    % default parameters
    schedule=1;
end

w=1;beta1=1;beta2=1;
dt1=0.001;
dt2=0.001;
optimum=0;

if schedule==1
    load("exponential_cooling_schedule.mat")
elseif schedule==2
    load("polynomial_cooling_schedule.mat")
end

gbest=1e10;
pbestx=zeros(n,N*p);
pbest=1e10*ones(1,N);
gbestx=zeros(n,N*p);

alpha_not_too_small=1;
all_iter=length(alpha_all);
tic
%% generate the centers randomly, and assign weights
[W2,pbest,pbestx,gbest,gbestx]=generate_centers_and_assign_w_initial_pbest_gbest(X,n,p,N,pbest,pbestx,gbest,gbestx);
for t=1:all_iter
    alpha=alpha_all(t);   
    it=0;
    count_gli=0;
    v=2*rand(n,N*p)-1;

    %% Due to computational accuracy, if alpha is annealed to a very small value, FCM degenerates to KM
    if alpha_not_too_small
        while true
            %% Multiple FCM modules for alternative clustering and determining the individual best.
            for nn_count=1:N
                W=W2(:,(nn_count-1)*p+1:nn_count*p);
                W=W';
                pre_obj=1e10;
                norm_flg=1;

                %% weight and theta iterations in FCM
                for iteration=1:600

                    mf = W.^alpha;                               % MF matrix after exponential modification
                    theta = mf*X./(sum(mf,2)*ones(1,size(X,2))); %new center

                    if (sum(isnan(theta))~=0)
                        alpha_not_too_small=0;
                        disp('NA');
                        norm_flg=0;
                        [W,obj]=calculate_k_means_obj_func_value(X,n,m,p,W);
                        break
                    end

                    dist = distfcm(theta, X);                  % fill the distance matrix
                    obj = sum(sum((dist.^2).*mf));
                    tmp = (dist+0.00001).^(-2/(alpha-1));      % calculate new U, suppose expo != 1
                    W = tmp./(ones(p, 1)*sum(tmp));

                    if pre_obj - obj < 1e-10
                        break
                    else
                        pre_obj=obj;
                    end
                end

                if norm_flg
                    [W,obj]=calculate_k_means_obj_func_value(X,n,m,p,W);
                end

                W2(:,(nn_count-1)*p+1:nn_count*p)=W;

                %% determine the individual best solution
                if obj<pbest(1,nn_count)
                    pbest(1,nn_count)=obj;
                    pbestx(:,(nn_count-1)*p+1:nn_count*p)=W;
                end
            end

            %% determine the group best solution
            [~,pc]=min(pbest);
            if pbest(1,pc)<gbest
                gbest=pbest(1,pc);
                gbestx=repmat(pbestx(:,(pc-1)*p+1:pc*p),1,N);
                count_gli=0;
                fprintf('objective function is updated: %0.8f\n',gbest)
                time_consume=toc;
                %disp(gbest)
            else
                count_gli=count_gli+1;
            end

            %% Termination Criterion - CNO
            if count_gli>M
                break
            end

            %% Compute the diversity and perform mutation if necessary
            Div=norm(W2-gbestx,'fro');
            Div=Div/(N*n*p);
            it=it+1;
            if Div < dt1
                disp('mutation!')
                a=exp(it/1000);%1-e
                psi=-2.5*a+5*a*rand(n,N*p);%-2.5a-2.5a
                et=(1/sqrt(a))*exp((-(psi/a).^2)/2).*cos((5*(psi/a)));
                larger_zero=et>0;
                lower_zero=et<0;
                W2=W2+larger_zero.*et.*(1-W2);
                W2=W2+lower_zero.*et.*W2;
                v=2*rand(n,N*p)-1;
                pbest=1e10*ones(1,N);
            else
                %% Reposition the initial weights by using a PSO rule
                v=w*v+beta1*rand(n,N*p).*(pbestx-W2)+beta2*rand(n,N*p).*(gbestx-W2);
                W2=W2+v;
                W2=min(1,max(0,W2));
            end
        end
    else
        %% alpha is annealed to a very small value, FCM degenerates to KM
        while true
            %% Multiple FCM modules for alternative clustering and determining the individual best.
            for nn_count=1:N
                W=W2(:,(nn_count-1)*p+1:nn_count*p);
                W=W';
                theta = W*X./(sum(W,2)*ones(1,size(X,2))); %new center

                index=isnan(sum(theta,2));
                if sum(index)~=0
                    theta(index,:)=X(randperm(n,sum(index)),:);
                end

                md = pdist2(X,theta,'squaredeuclidean');
                [md_min,label]=min(md,[],2);
                pre_obj=sum(md_min);

                while 1
                    for l = 1:p
                        theta(l, :) = mean(X(label == l, :),1);
                    end
                    md = pdist2(X,theta,'squaredeuclidean');
                    [md_min,label]=min(md,[],2);
                    obj=sum(md_min);

                    if pre_obj-obj<1e-10
                        break
                    else
                        pre_obj=obj;
                    end
                end
                W=zeros(n,p);
                for i=1:n
                    W(i,label(i))=1;
                end

                W2(:,(nn_count-1)*p+1:nn_count*p)=W;
                %% determine the individual best solution
                if obj<pbest(1,nn_count)
                    pbest(1,nn_count)=obj;
                    pbestx(:,(nn_count-1)*p+1:nn_count*p)=W;
                end
            end

            %% determine the group best solution
            [~,pc]=min(pbest);
            if pbest(1,pc)<gbest
                gbest=pbest(1,pc);
                gbestx=repmat(pbestx(:,(pc-1)*p+1:pc*p),1,N);
                count_gli=0;
                fprintf('objective function is updated: %0.8f\n',gbest)
                time_consume=toc;
                %disp(gbest)
            else
                count_gli=count_gli+1;
            end

            %% Termination Criterion - CNO
            if count_gli>M
                break
            end

            %% Compute the diversity and perform mutation if necessary
            Div=norm(W2-gbestx,'fro');
            Div=Div/(N*n*p);
            it=it+1;
            if Div < dt2
                disp('mutation!')
                a=exp(it/1000);%1-e
                psi=-2.5*a+5*a*rand(n,N*p);%-2.5a-2.5a
                et=(1/sqrt(a))*exp((-(psi/a).^2)/2).*cos((5*(psi/a)));
                larger_zero=et>0;
                lower_zero=et<0;
                W2=W2+larger_zero.*et.*(1-W2);
                W2=W2+lower_zero.*et.*W2;
                v=2*rand(n,N*p)-1;
                pbest=1e10*ones(1,N);
            else
                %% Reposition the initial weights by using a PSO rule
                v=w*v+beta1*rand(n,N*p).*(pbestx-W2)+beta2*rand(n,N*p).*(gbestx-W2);
                W2=W2+v;
                W2=min(1,max(0,W2));
            end
        end
    end
    [W2,pbest,pbestx,gbest,gbestx]=generate_centers_and_assign_w_initial_pbest_gbest(X,n,p,N,pbest,pbestx,gbest,gbestx);
    if gbest<=optimum
        break
    end
end
all_time=toc;
disp('end')

W=gbestx(:,1:p);
[~,label]=max(W,[],2);
end

function [W,obj]=calculate_k_means_obj_func_value(X,n,m,p,W)
[~,label]=max(W);
theta=zeros(p,m);
for l = 1:p
    theta(l, :) = mean(X(label == l, :),1);
end
index=isnan(sum(theta,2));
if sum(index)~=0
    theta(index,:)=X(randperm(n,sum(index)),:);
end

md = pdist2(X,theta,'squaredeuclidean');
[md_min,label]=min(md,[],2);
pre_obj=sum(md_min);

while 1
    for l = 1:p
        theta(l, :) = mean(X(label == l, :),1);
    end
    md = pdist2(X,theta,'squaredeuclidean');
    [md_min,label]=min(md,[],2);
    obj=sum(md_min);

    if pre_obj-obj<1e-10
        break
    else
        %disp(pre_obj-obj)
        pre_obj=obj;
    end
end
W=zeros(n,p);
for i=1:n
    W(i,label(i))=1;
end
end

function [W2,pbest,pbestx,gbest,gbestx]=generate_centers_and_assign_w_initial_pbest_gbest(X,n,p,N,pbest,pbestx,gbest,gbestx)
sz=[n,p];
row=[1:n];

for nn_count=1:N
    sam=randperm(n,p);
    theta=X(sam,:);

    md = pdist2(X,theta,'squaredeuclidean');
    [~,label]=min(md,[],2);
    [obj,~,label] = label2wd_2(X,label,p);
    pre_obj=obj;
    while 1
        for l = 1:p
            theta(l, :) = mean(X(label == l, :),1);
        end
        md = pdist2(X,theta,'squaredeuclidean');
        [md_min,label]=min(md,[],2);
        obj=sum(md_min);
        if pre_obj-obj<1e-10
            break
        else
            %disp(pre_obj-obj)
            pre_obj=obj;
        end
    end

    W=zeros(n,p);
    W(sub2ind(sz,row,label'))=1;
    pbest(1,nn_count)=obj;
    pbestx(:,(nn_count-1)*p+1:nn_count*p)=W;

end
W2=pbestx;

[best,pc]=min(pbest);
if best<gbest
    gbest=best;
    gbestx=repmat(pbestx(:,(pc-1)*p+1:pc*p),1,N);
    fprintf('objective function is updated in initialization: %0.8f\n',gbest)
end

for i=1:size(W2,2)
    if all((W2(:,i))==0)
        W2(randperm(n,p),i)=1;
    end
end
end