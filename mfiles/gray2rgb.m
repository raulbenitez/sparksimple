function [nI] = gray2rgb(I,map)



nI=zeros(size(I,1),size(I,2),3);
sc0=zeros(size(I,1),size(I,2));

rg=[min(min(I)):max(max(I))];
SM=size(map,1);
for ii=1:SM
    try
        if(ii==SM)
    qui=find(I>=rg(ii));        
        else
    qui=find(I==rg(ii));
        end
    sc=sc0;
    
    sc(qui)=map(ii,1);
    nI(:,:,1)=nI(:,:,1)+sc;
    
    sc(qui)=map(ii,2);
    nI(:,:,2)=nI(:,:,2)+sc;
    
    sc(qui)=map(ii,3);
    nI(:,:,3)=nI(:,:,3)+sc;
    end
end



end