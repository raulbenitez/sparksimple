function [mask] = cellMask4(mk)

mk0=mk;
mk=double(mk);
mk=(mk-min(min(mk)))/(max(max(mk))-min(min(mk)));
if((max(max(mk0))>256)||(max(max(mk0))<=1))
ed=[0:.001:1];
else
ed=[0:.004:1];
end
h=histc(reshape(mk,1,numel(mk)),ed);
wind=ceil(length(ed)/16);
h=filtfilt(ones(1,wind)/wind,1,h);

%figure;plot(ed,h)

hd=diff(h);
P=hd(1:end-1).*hd(2:end);
zr=find(P<=0)+1;
hd2=diff(hd);

% aux=zr(find(hd2(zr)>=0));
% figure,plot(ed,h);hold on;plot(ed(aux),h(aux),'r+');
%figure;plot(h),hold on,plot(P);plot(hd2)
thresh=[];
% try
% thresh=[thresh ed(zr(find(hd2(zr)>=0,1,'first')))];
% end
try
thresh=[thresh ed(zr(find(hd2(zr+1)>=0,1,'first')))];
end
try
    aux=reshape(mk,numel(mk),1);
    aux=sort(aux);
IDX=kmeans(aux,2);
thresh=[thresh mean(sort([mean(aux(IDX==1)) mean(aux(IDX==2))]))];
end
if(isempty(thresh)),thresh=.5;else thresh=mean(thresh);end

mask=mk>thresh;

% keep largest
BLUE=mask;
BLUE=bwlabel(BLUE);
kk=regionprops(BLUE,'area');
kk=cell2mat(struct2cell(kk));
[v,ik]=max(kk);
BLUE(BLUE~=ik)=0;
BLUE(BLUE~=0)=1;

%remove holes
L=bwlabel(1-BLUE);
kk=regionprops(L,'PixelList');
for ixi=1:length(kk)
   vals=kk(ixi).PixelList; 
   if(isempty([find(vals(:,1)==1); find(vals(:,2)==1);find(vals(:,2)==size(L,1));find(vals(:,1)==size(L,2))]))
       BLUE(L==ixi)=1;
   end
end

mask=BLUE;

end