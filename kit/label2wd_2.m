function [wd,theta,label] = label2wd_2(data,label,p)
m=size(data,2);
theta=zeros(p,m);
for l = 1:p
    theta(l, :) = mean(data(label == l, :));
end

md = pdist2(data,theta,'squaredeuclidean');
[md_min,label]=min(md,[],2);            
wd=sum(md_min);   
                 
end
