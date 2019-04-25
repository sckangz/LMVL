function [result]=multiviewsemi(X,y,r,labelind,beta,gamma)%y用于衡量精度
mv=size(X,2);
[m,num]=size(X{1}');
H=rand(r,num);
c=length(unique(y));%c为簇种类
S=eye(num);%S的初始化
k=10;
% initialize F
D = diag(sum(S));
L = D - S;
  ulabel= setdiff(1:num,labelind);%index of unlabeled data

    ll=length(labelind);
    Fl=zeros(ll,c);
    for i=1:ll
        Fl(i,y(labelind(i)))=1;
    end
        Fu=zeros(num-ll,c);
        rr = zeros(num,1);
   %     alpha=100;
        F=[Fl;Fu];
%     for i = 1:num
%     di = beta*distH1(i,2:k+2)+gamma*distF1(i,2:k+2);
% %     rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
%    id = idh(i,2:k+2);
%      S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
%     end
options = optimset( 'Algorithm','interior-point-convex','Display','off');  for j=1:50
    for i=1:mv
        A=X{i}';
for ij=1:m
    Wv(ij,:) = quadprog(2*H*H',-2*A(ij,:)*H',[],[],[],[],zeros(num,1),[],[],options); 
end
% size(Wv)
wv{i}=Wv;      
        
        M(:,:,i)=wv{i}'*wv{i};
        N(:,:,i)=(wv{i})'*X{i}';
    end
    H=sylvester(sum(M,3),beta*L,sum(N,3));%sum中dim是否使用正确？？？
    F_old=F;
    dist_f=L2_distance_1(F',F');
    dist_h=L2_distance_1(H,H);
    [distH1,idh]=sort(dist_h,2);
    [distF1,idf]=sort(dist_f,2);
    
    for i = 1:num
    di = beta*distH1(i,2:k+2)+gamma*distF1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
     end

alpha = mean(rr)/(2*beta);

    for i=1:num
       idha0=idh(i,2:k+1);
       dfi=dist_f(i,idha0);
       dhi=dist_h(i,idha0);
       ad=-(beta*dhi+gamma*dfi)/(4*alpha*beta);
   S(i,idha0)=EProjSimplex_new(ad);
   
   end
S=S-diag(diag(S));

    S = (S+S')/2;                                                        %update F
    D = diag(sum(S));
    L = D-S;
   uu=zeros(num-ll,num-ll);
        for ii=1:(num-ll)
            for jj=1:(num-ll)
                uu(ii,jj)=L(ulabel(ii),ulabel(jj));
            end
        end
        ul=zeros(num-ll,ll);
        for ii=1:(num-ll)
            for jj=1:ll
                ul(ii,jj)=L(ulabel(ii),labelind(jj));
            end
        end
        Fu=-uu\(ul*Fl);
            F=[Fl;Fu];
  if ((j>1)&&(norm(F-F_old,'fro')<norm(F,'fro')*1e-5))
        break
    end
%     evs(:,j+1) = ev;
%         %update lambda
%     thre = 1*10^-8;
%     fn1 = sum(ev(1:c));                                                %update lambda
%     fn2 = sum(ev(1:c+1));
%     if fn1 > thre
%         gamma = 2*gamma;
%     elseif fn2 < thre
%         gamma = gamma/2;  F = F_old;
%     else
%         break;
%     end
end
[ur,uc]=size(ulabel);
    [max_value,max_ind] = max(Fu,[],2);
    cnt = 0;
    for i = 1:uc
        if max_ind(i) == y(ulabel(i))
            cnt = cnt+1;
        end
    end
    result = cnt/uc;
   
% L=(S+S')/2;
% actual_ids = spectral_clustering(L, c);
% [result] = ClusteringMeasure( actual_ids,y);
% [clusternum, predy]=graphconncomp(sparse(S)); predy = predy';
% if clusternum ~= c
%     sprintf('Can not find the correct cluster number: %d', clusternum)
% end;
% result = ClusteringMeasure(y, predy);
    
