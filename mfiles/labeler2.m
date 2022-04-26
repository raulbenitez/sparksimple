function [im2] = labeler2(im,msk)

[a,b]=size(im);
safe=6;% radi del q creiem ocuparan (per no solaparse)
labs=unique(im);labs(labs==0)=[];

pp=zeros(length(labs),4);
dd=bwdist(1-msk);
z0=zeros(size(im));
im2=z0;imsf=z0;
for ii=1:length(labs)

    
    pp(ii,1)=labs(ii);
    [y2,x2]=find(im==labs(ii));
if(~isempty(y2))
   pp(ii,2)=y2(1);pp(ii,3)=x2(1); 
   radi=ceil(dd(y2(1),x2(1)));
   zz=z0;
   zz(y2(1),x2(1))=1;
   zz=imdilate(zz,strel('disk',double(radi)));
   reg=2*zz-msk==2;% figure,imagesc(2*zz-msk)
          L=bwlabel(reg);aux=unique(L);aux(aux==0)=[];
          if(isempty(aux))
              cand=[];
          else
          salva=ceil(rand(1)*length(aux));if(salva==0),salva=1;end
          L(L~=aux(salva))=0;
       cand=find(L>0);% figure,imagesc(reg)
          end
   while(isempty(cand))
       zz=imdilate(zz,[1 1 1;1 1 1;1 1 1]);
       reg=2*zz-msk==2;
  L=bwlabel(reg);aux=unique(L);aux(aux==0)=[];
          if(isempty(aux))
              cand=[];
          else
          salva=ceil(rand(1)*length(aux));if(salva==0),salva=1;end
          L(L~=aux(salva))=0;
       cand=find(L>0);% figure,imagesc(reg)
          end
   end
   [y,x]=ind2sub([a,b],cand);
   x=round(mean(x));y=round(mean(y));
   pos=[y2(1),x2(1)];
   v=[y,x]-pos;
  % vn=v/norm(v);
   % quants cops el vector director fins a interseccio limits
   %%% top    left      bot   right
    n=[([1,1]-pos)./v,([a,b]-pos)./v];
   n(n<0)=max([a,b]);
   [~,q]=min(n);
   add=[0,0];
   if(q==1),add=[safe+1,0];end
   if(q==2),add=[0,safe+1];end
   if(q==3),add=[-safe-1,0];end
   if(q==4),add=[0,-safe-1];end  
   posf=pos+round(v*n(q))+add;
   imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)=imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)+1;
   % figure,imagesc(imsf),axis image
   inters=find(imsf>1);
   if(~isempty(inters))
       try
       imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)=imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)-1;
       [y,x]=ind2sub([a,b],inters);
       vv=posf-mean([y,x]);N=norm(vv);
       if(N==0),vv=add*2;else,vv=vv/N;vv=ceil([max(y)-min(y)+1,max(x)-min(x)+1]).*vv;end       
       posf=posf+round(vv);
       imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)=imsf(posf(1)-safe:posf(1)+safe,posf(2)-safe:posf(2)+safe)+1;
       end
   end
   %%%%%
    [lesx,lesy]=puntsEnmig([pos(2),posf(2)],[pos(1),posf(1)]);
    q=find(lesx<1);lesy(q)=[];lesx(q)=[];
    q=find(lesy<1);lesx(q)=[];lesy(q)=[];
    for jj=1:length(lesy)
         im2(lesy(jj),lesx(jj))=1;
    end
    hl='center';vl='middle';
 
   im2=textIm(posf(2),posf(1),num2str(labs(ii)),im2,'horizontalalignment',hl,'verticalalignment',vl,'blending','off','background','box','bgcolor',0,'bdcolor',1);
% figure,imagesc(im2),axis image
end
end
% [im2]=labelsEquidist(pp,im);
% figure(9),imagesc(im2+msk),axis image
% 
% try,
%     [im2]=labelsRadial(pp,im);
% catch  
% try
%     [im2]=labelsStraight(pp,im);
% catch
% [im2] = labelsEquidist(pp,im);
% end
% end

end