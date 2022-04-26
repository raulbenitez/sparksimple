function [] = batchFishResults(fol)


tagRyR='c1';
tagSpk='c2';
R=[];N=[];Z=[];DATA=[];%P=[];
FAN=zeros(57,57);% atenció, 57 pel fet q són tots els experiments amb el 
% mateix DX i podem composar fancies, sino no es pot fer

S=dir(fol);

for jj=3:length(S)
    
 SR=[dir([fol '/' S(jj).name '/rgb/*' tagRyR '*.jpg']); dir([fol '/' S(jj).name '/rgb/*' tagRyR '*.tif'])];
 SS=[dir([fol '/' S(jj).name '/rgb/*' tagSpk '*.jpg']); dir([fol '/' S(jj).name '/rgb/*' tagSpk '*.tif'])];
 
 
 if(length(SR)~=length(SS))
     error('Number of images in each channel is not the same!');
 else
 
     
     
     
 Rfol=[fol '/' S(jj).name '/rgb/Results'];
 fol2=[fol '/' S(jj).name '/rgb/'];
 %Rfol2=[fol '/sarco'];
 
 ensenya(['Folder ' fol2]);
  % process RyR individually

    load([Rfol '/zData2.mat']);
   
   dades=guardaDades(spkF,fol2,Rfol);
    DATA=[DATA;dades];
 %fancy
 FAN=FAN+fancy;
 
 
 z=zeros(length(SR),1);
 res=getRes(S(jj).name);
  for kk=1:length(SR)
  load([Rfol '/' SR(kk).name(1:end-4) '/zData.mat']);
   z(kk)=estimateZlineDist(pos,res);%close all;
   R=[R;rads/res];
   N=[N;nn'/res];
  % P=[P;pdist(pos*res)'];
  end
   z(z==0)=[];
   Z=[Z;z];
 end

end

Z=Z/res;
qq=[.5];
radis=quantile(R,qq);
t = linspace(0,2*pi); 

% [ka,kb]=find(FAN>0);
% desp=min([min(ka) size(FAN,1)-max(ka)+1 min(kb) size(FAN,2)-max(kb)+1]);
% nr=desp:size(FAN,1)-desp+1;
% FAN=FAN(desp:size(FAN,1)-desp+1,desp:size(FAN,1)-desp+1);
% outs=[find(ypos<nr(1)) find(ypos>nr(end))];
% ypos(outs)=[];yvals(outs)=[];
% ypos=ypos-nr(1)+1;
 [ka,kb]=size(FAN);

figure;set(gcf,'position',[1326         678         683         649],'color','w');
aux0=110;
CM=jet(aux0);CM(1,:)=[0 0 0];
aux=round(aux0*.9);
CM(aux+1:256,:)=[(linspace(CM(aux,1),4,256-aux))' linspace(CM(aux,2),1,256-aux)' linspace(CM(aux,3),2,256-aux)'];
CM(CM>1)=1;
imagesc(FAN);set(gca,'fontsize',14);colormap(CM);hold on;
for ii=1:length(radis)
    if(ceil(length(radis)/2)==ii)
       sz=2; 
    else
     sz=1;   
    end
    h=plot(ceil(ka/2)+radis(ii)*cos(t),ceil(ka/2)+radis(ii)*sin(t),'y');set(h,'linewidth',sz);
   % h=text(ceil(ka/2)+radis(ii),ceil(ka/2),[' ' num2str(round(qq(ii)*100)) '%']);set(h,'color','w','horizontalalignment','left','fontsize',8);
    if(ii==1)
        %est=(.26/getRes(S(3).name));
            %h=plot(ceil(ka/2)+est*cos(t),ceil(ka/2)+est*sin(t),'r');%set(h,'linewidth',2);
            est=mean(N);
            h=plot(ceil(ka/2)+est*cos(t),ceil(ka/2)+est*sin(t),':m');
            h=plot(ceil(ka/2)+mean(Z)*cos(t),ceil(ka/2)+mean(Z)*sin(t),':c');
%[h,ha]=legend('RyR mean apparent size','RyR estimated size','RyR nearest neighbour','Z-line distance');
[h,ha]=legend('RyR mean apparent size','RyR nearest neighbour','Z-line distance');
set(h,'color','k');
for jj=1:3
   set(ha(jj),'color','w'); 
end
    end
end
%est=(.26/getRes(S(3).name));
%h=plot(ceil(ka/2)+est*cos(t),ceil(ka/2)+est*sin(t),'r');set(h,'linewidth',2);
set(gca,'xtick',ypos,'xticklabel',yvals,'ytick',ypos,'yticklabel',yvals);
title(['RyR-spark colocalization (N=' num2str(size(DATA,1)) ')'],'fontsize',14);xlabel('distance to RyR (um)');axis equal;colorbar;

% grid lines
for ii=1:length(ypos)
   h=line([0.5 ka+.5],[ypos(ii) ypos(ii)]);set(h,'color',[.5 .5 .5],'linestyle',':');
   h=line([ypos(ii) ypos(ii)],[0.5 ka+.5]);set(h,'color',[.5 .5 .5],'linestyle',':');
end

[a,b]=find(FAN>0);
a1=[];b1=[];off=round(size(FAN,1)/2);
for ii=1:length(a)
   quants=FAN(a(ii),b(ii));
   for jj=1:quants
      a1=[a1 a(ii)-off];
      b1=[b1 b(ii)-off];
   end
end


dist=sqrt((a1.*a1)+(b1.*b1));
dist=sort(dist);
ratios=[.5 .75 .9 .99];
distancies=zeros(size(ratios));
for ii=1:length(ratios)
posicio=round(ratios(ii)*length(dist));
distancies(ii)=round(100*dist(posicio)*getRes(S(3).name))/100;

h=text(1,2*ii,[num2str(round(ratios(ii)*100)) '% @ d < ' num2str(distancies(ii)) 'um']);
set(h,'color','w','fontsize',14);
end





%saveWysiwyg(gcf,['C:\Users\Austin Powers\Desktop\SPK3.png']);
%close all;
pintaDades(DATA);
ensenya(['Surviving sparks: ' num2str(size(DATA,1)) '.'])

end