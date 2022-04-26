function L2 = proximityFilter(IX,dm,tm,DX,DT)

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

reg=cell(elems,1);
for ii=1:elems
   reg{ii}=ii; 
end
for ii=1:length(outs)
    
    ind=[ind1(outs(ii)) ind2(outs(ii))];
    
    reg{ind(1)}=[reg{ind(1)} ind(2)];
     reg{ind(2)}=[reg{ind(2)} ind(1)];
end
for ii=1:elems
   reg{ii}=sort(unique(reg{ii})); 
end


L0=L;
L2=L;L2(L2>0)=0;
%cont=elems;
for ii=1:elems
reps=reg{ii}; 
vals=reps;
for jj=1:length(reps)
    aux=IX(L==reps(jj));
    vals(jj)=aux(1);
    %L(L==reps(jj))=0;
end
aux=round(mean(c02(reps,:),1));
L2(aux(1),aux(2),aux(3))=round(mean(vals));
    
    
end








end