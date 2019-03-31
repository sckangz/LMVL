load bbc_seg14of4.mat
warning off
r=[ 20 50] ;
beta=[ 0.01 0.1 1 ];
gamma=[ 0.1 1 10   ];
sonarresult=zeros(length(r),1);
rate=0.5 ;%rate of labeled data
% for i=1:size(data,1)
% dist = max(max(data{i})) - min(min(data{i}));
% m01 = (data{i} - min(min(data{i})))/dist;
% data{i} = 2 * m01 - 1;
% end
num=size(data{1},1);
c=length(unique(labels)); % number of class
numperc=floor(num/c); % number of data per class
labelperc=floor(rate*numperc); % number of labeled data per class
fid=fopen('latentgraphresults.txt','a');
fprintf(fid,'%25s','bbc_seg14of4.mat');
fprintf(fid,'%12s %12.6f\r\n','rate is',rate);
fprintf(fid,'%12s %12s %12s %12s\r\n','r','beta','gamma','result');
for t=1:20
labelindperc=sort(randperm(numperc,labelperc)); % index of labeled data selected
labelind=[]; % labelind: index of known label
for i=1:c
    labelind=[labelind labelindperc+(i-1)*numperc];
end
lres=0;
for i=1:length(r)
     for j=1:length(beta)
         for k=1:length(gamma) 
          fprintf('params%12.6f%12.6f%12.6f\n',r(i),beta(j),gamma(k));
           [result]=multiviewsemi(data,labels,r(i),labelind,beta(j),gamma(k))
            fprintf(fid,'%12.6f %12.6f %12.6f %12.6f\r\n',r(i),beta(j),gamma(k),result);
            if (result>lres)
            lres=result;
            end
         end
     end
end
 acc(t)=lres;
end


 mean(acc)
 std(acc)
 fprintf(fid,'%12s %12.6f\r\n','mean',mean(acc) );
 fprintf(fid,'%12s %12.6f\r\n','std',std(acc) );
 fclose(fid);