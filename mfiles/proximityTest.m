function [] = proximityTest(IX,dm,tm,DX,DT)

[a0 b0 c0]=size(IX);
IXX=double(IX>0);
L=bwlabeln(IXX);
elems=max(max(max(L)));
c2=zeros(elems,3);

for ii=1:elems
        [a1,b1]=find(L==ii);
        if(length(a1)>1),a1=a1(1);b1=b1(1);end
     [b1,c1]=ind2sub([b0,c0],b1);
    c2(ii,1:3)=[a1 b1 c1];
    
end



c02=c2;
for ii=1:size(c2,1)
c2(ii,:)=c2(ii,:).*[DX DX DT];
end
dist=(pdist(c2(:,1:2)));
time=(pdist(c2(:,3)));
same=[];

ind1=[];ind2=[];
for ii=1:size(c2,1)-1
    for jj=ii+1:size(c2,1)
        ind1=[ind1 ii];
        ind2=[ind2 jj];
    end
end

outs=find((dist<dm)&(time<tm));
if(numel(outs)>0),
ensenya(['Warning: Found ' num2str(numel(outs)) ' distances out of ' num2str(numel(dist)) ' (' num2str(elems) ' elements).'],superjet(1,'a'));

end
end