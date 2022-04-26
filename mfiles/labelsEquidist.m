function [im2] = labelsEquidist(pp,im)

[a,b]=size(im);
labs=pp(:,1);
[aux,ord]=sort(pp(:,4),'descend');
pp2=pp(ord,:);
im2=zeros(size(im));
peri=[a b a b];
peri=round(length(labs)*cumsum(peri)/sum(peri));
%figure;imagesc(mask);colormap(gray(2));hold on;axis equal
inte=a/peri(1);ij=1;
for ii=1:peri(1)

%     plot(pp2(ii,3),pp2(ii,2),'+r');    
%     plot([pp2(ii,3),1],[pp2(ii,2),ij],'r');
    
    try im2(pp2(ii,2)-1:pp2(ii,2)+1,pp2(ii,3))=1;end
    try im2(pp2(ii,2),pp2(ii,3)-1:pp2(ii,3)+1)=1;end
    T=escriuFrase(num2str(pp2(ii,1)),12);T=double(T);T(15:end,:)=[];[sv,sh]=size(T);
   mgh=round(sh/2);mgv=round(sv/2);
   x2=1+mgh+1;y2=round(ij+mgv+1);
      try im2(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh-1)=1;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh+sh+1)=1;end
      try im2(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
      try im2(y2-mgv+sv+1,x2-mgh-1:x2-mgh+sh+1)=1;end
      if(ii==1),y2=y2-mgv+sv+1;end
      if(ii==peri(1)),y2=y2-mgv-1;end
   [lesx,lesy]=puntsEnmig([x2-mgh+sh+1,pp2(ii,3)],[y2,pp2(ii,2)]);
   for jj=1:length(lesy)
       im2(lesy(jj),lesx(jj))=1;
   end
   
    ij=ij+inte;
    
end
if(exist('mgh','var')==0),
ij=1;
else
ij=x2-mgh+sh+1+1;
end
inte=(b-2*ij)/(peri(2)-peri(1));
for ii=peri(1)+1:peri(2)

%     plot(pp2(ii,3),pp2(ii,2),'+g');    
%     plot([pp2(ii,3),ij],[pp2(ii,2),a],'g');
    try im2(pp2(ii,2)-1:pp2(ii,2)+1,pp2(ii,3))=1;end
    try im2(pp2(ii,2),pp2(ii,3)-1:pp2(ii,3)+1)=1;end
    T=escriuFrase(num2str(pp2(ii,1)),12);T=double(T);T(15:end,:)=[];[sv,sh]=size(T);
   mgh=round(sh/2);mgv=round(sv/2);
   x2=round(ij+mgh+1);y2=a-mgv-1;
      try im2(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh-1)=1;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh+sh+1)=1;end
      try im2(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
      try im2(y2-mgv+sv+1,x2-mgh-1:x2-mgh+sh+1)=1;end
   [lesx,lesy]=puntsEnmig([x2,pp2(ii,3)],[y2-mgv-1,pp2(ii,2)]);
   for jj=1:length(lesy)
       im2(lesy(jj),lesx(jj))=1;
   end


ij=ij+inte;
    
end
mgh=0;clear mgh;
inte=a/(peri(3)-peri(2));ij=a;
for ii=peri(2)+1:peri(3)

%     plot(pp2(ii,3),pp2(ii,2),'+b');    
%     plot([pp2(ii,3),b],[pp2(ii,2),ij],'b');

   try im2(pp2(ii,2)-1:pp2(ii,2)+1,pp2(ii,3))=1;end
    try im2(pp2(ii,2),pp2(ii,3)-1:pp2(ii,3)+1)=1;end
    T=escriuFrase(num2str(pp2(ii,1)),12);T=double(T);T(15:end,:)=[];[sv,sh]=size(T);
   mgh=round(sh/2);mgv=round(sv/2);
   x2=b-mgh-1;y2=round(ij-mgv-1);
      try im2(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh-1)=1;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh+sh+1)=1;end
      try im2(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
      try im2(y2-mgv+sv+1,x2-mgh-1:x2-mgh+sh+1)=1;end
            if(ii==peri(3)),y2=y2-mgv+sv+1;end
      if(ii==peri(2)+1),y2=y2-mgv-1;end
   [lesx,lesy]=puntsEnmig([x2-mgh-1,pp2(ii,3)],[y2,pp2(ii,2)]);
   for jj=1:length(lesy)
       im2(lesy(jj),lesx(jj))=1;
   end

ij=ij-inte;
    
end
if(exist('mgh','var')==0),
ij=b;
else
ij=x2-mgh-1-1;
end

inte=(b-2*(b-ij))/(length(labs)-peri(3));
for ii=peri(3)+1:length(labs)

%     plot(pp2(ii,3),pp2(ii,2),'+m');    plot([pp2(ii,3),ij],[pp2(ii,2),1],'m');

    try im2(pp2(ii,2)-1:pp2(ii,2)+1,pp2(ii,3))=1;end
    try im2(pp2(ii,2),pp2(ii,3)-1:pp2(ii,3)+1)=1;end
    T=escriuFrase(num2str(pp2(ii,1)),12);T=double(T);T(15:end,:)=[];[sv,sh]=size(T);
   mgh=round(sh/2);mgv=round(sv/2);
   x2=round(ij-mgh-1);y2=1+mgv+1;
      try im2(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh-1)=1;end
      try im2(y2-mgv:y2-mgv+sv,x2-mgh+sh+1)=1;end
      try im2(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
      try im2(y2-mgv+sv+1,x2-mgh-1:x2-mgh+sh+1)=1;end
   [lesx,lesy]=puntsEnmig([x2,pp2(ii,3)],[y2+mgv+1,pp2(ii,2)]);
   for jj=1:length(lesy)
       im2(lesy(jj),lesx(jj))=1;
   end

    ij=ij-inte;
    
end


end