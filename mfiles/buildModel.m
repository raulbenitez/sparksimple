function [] = buildModel(folder,resultsFolder)

folder(folder=='\')='/';
if(folder(end)~='/'),folder=[folder '/'];end
resultsFolder(resultsFolder=='\')='/';
if(resultsFolder(end)~='/'),resultsFolder=[resultsFolder '/'];end
if(exist(resultsFolder,'dir')==0),mkdir(resultsFolder);end
 endFolder=folder(find(folder(1:end-1)=='/',1,'last')+1:end-1);
    resultsFolder=[resultsFolder endFolder];
    
    
  
    load([resultsFolder '/3Ddata.mat']);
    
    resolutions=getRes(folder);%microns per pixel
resolution=resolutions(1);
try
    zdim=resolutions(2);%microns per pixel
end
    
Q(:,1:2)=Q(:,1:2)*resolution;
Q(:,3)=Q(:,3)*zdim;
    
   
    
    
    dis=pdist(Q(:,1:3));
    dis=squareform(dis);
        dis(dis==0)=9999;
        
        [d o]=sort(dis,2);

N=2;



figure(33);hold on;
    
    for ii=1:size(Q,1)
        if((Q(ii,3)<26*zdim)&&(Q(ii,3)>14*zdim))
        plot3(Q(ii,2),Q(ii,1),Q(ii,3),'ro');
        for jj=1:N
            if((Q(o(o(ii,jj),jj),3)<40*zdim)&&(Q(o(ii,jj),3)>14*zdim))
           h=line([Q(ii,2) Q(o(ii,jj),2)],[Q(ii,1) Q(o(ii,jj),1)],[Q(ii,3) Q(o(ii,jj),3)]); 
            end
        end
        end
    end
    axis equal;
%     
% cc=0;
%      for ii=1:size(Q,1)
%      if((Q(ii,3)<40*zdim)&&(Q(ii,3)>14*zdim))
%          cc=cc+1;
%     v(cc,1:3)=[Q(ii,2)-Q(o(ii,jj),2),Q(ii,1)-Q(o(ii,jj),1),Q(ii,3)-Q(o(ii,jj),3)];
%      end
%      end
%      VT=[0 0 0];
%      for ii=1:cc
%          for jj=1:cc
%              VT=VT+(v(cc,:)/norm(v(cc,:)));
%          end
%      end
%     
%     


end
