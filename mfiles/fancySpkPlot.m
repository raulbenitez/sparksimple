function [] = fancySpkPlot(sigspk,catspk,DT)




pos=[129         103        1835         364]; % figure position
szx=.9;szx0=.05;% horizontal axes size and offset both panels
szy1=.7;szy01=0.25;% vertical axes size and offset top panel
szy2=.15;szy02=0.02;% vertical axes size and offset bot panel
rt=0;% remove text

fs=22;% font size
colorl1=[.5 .5 .5];% vertical lines 
colorl2=[.5 .5 .5];% tilted lines 
invertPatch=1;% shade relevant or irrelevant area
colorf0=[.9 .9 .9]; % figure background color
colora0=[1 1 1]; % axes background color
colorp1=[.90 .90 .90]; % top panel patch
colorp2=[.90 .90 .90]; % between panels patch
colorl3=[0 0 0];%patch border lines
lins=0;% lines in top panel
nums=0;% numbers in top panel
linorpatch=2;% 0 for lines 2 for patch 1 for both joining axes
neither=0;% neither lines nor patch (ignores previous parameter)
twolines=0;% draw first and last lines regardless of previous parameter
centersborders=1;% 1 for centers 2 for borders
patchborder=1;% patch iner border line
cm=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];cm(1,:)=[0 0 0];
 
or=find(sigspk==max(sigspk))-5;
fi=find(sigspk==max(sigspk))+12;
%=round((size(catspk,2)+2)/(size(catspk,1)+2));
ti=1:length(sigspk);
ts=ti*DT;
figure;clf;set(gcf,'color',colorf0,'position',pos);
xl=[0.5*DT,(length(sigspk)+0.5)*DT];
yl(1)=floor(min(sigspk));
yl(2)=ceil(max(sigspk));
    axes('position',[szx0 szy02+szy2 szx szy01-szy02-szy2],'xtick',[],'ytick',[]);hold on;
    if(invertPatch==0)
        set(gca,'color',colorf0);
    else
         set(gca,'color',colora0);
    end
     xlim(xl);ylim(yl);
      if((linorpatch>=1)&&(neither==0))
    if(invertPatch==0)
    patch([xl(1)/DT xl(2)/DT fi+.5 or-.5]*DT,[yl(1) yl(1) yl(2) yl(2)],colorp2,'edgecolor','none');
    else
    patch([xl(1)/DT or-.5 xl(1)/DT]*DT,[yl(1) yl(2) yl(2)],colorp2,'edgecolor','none');
    patch([ xl(2)/DT xl(2)/DT fi+.5]*DT,[ yl(1) yl(2) yl(2)],colorp2,'edgecolor','none');
    end          
      end
    HM=axes('position',[szx0 szy01 szx szy1]);hold on;
    set(HM,'color',colora0);
    if((neither==0))
    if((invertPatch==0))
    patch([or-.5 fi+.5 fi+.5 or-.5]*DT,[yl(1) yl(1) yl(2) yl(2)],colorp1,'edgecolor','none');
    else
    patch([xl(1)/DT or-.5 or-.5 xl(1)/DT]*DT,[yl(1) yl(1) yl(2) yl(2)],colorp1,'edgecolor','none','FaceVertexAlphaData',.5);
    patch([fi+.5 xl(2)/DT xl(2)/DT fi+.5]*DT,[yl(1) yl(1) yl(2) yl(2)],colorp1,'edgecolor','none','FaceVertexAlphaData',.1);
    end
    end
      for ixi=1:length(ti)
          ix=ixi*DT; 
          if(lins==1)
           h=line([ix,ix],[2*yl(1),sigspk(ixi)]);set(h,'color',colorl1,'linewidth',1);
          end
          if(nums==1)
           h=text([ix],[sigspk(ixi)],[num2str(ixi)]);set(h,'color',colorl1,'linewidth',1,'horizontalalignment','center','verticalalignment','bot');
          end
          %%%% h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-.5)*szx/length(ti)],[szy01+(szy1)*(ts(ixi)-yl(1))/diff(yl) szy02+szy2]);set(h,'color',[.2 .2 .2],'linewidth',.5);
      end
      if((linorpatch<=1)&&(neither==0))
      for ixi=or:fi
          ix=ixi*DT;
          if(centersborders==1)
            h=annotation('line',[szx0+(ixi)*szx/length(ti) szx0+(ixi-or+.97)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl2,'linewidth',.5);  
          end
          if(centersborders==2)
           h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-or+.5)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl2,'linewidth',.5);    
          end
      end
      end
      if(twolines==1),
          ixi=or;ix=ixi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-or+.5)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl2,'linewidth',.5);
          ixi=fi;ix=fi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-or+.5)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl2,'linewidth',.5);
      end
      if((patchborder==1)&&(neither==0)),
          ixi=or-.5;ix=ixi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-or+.5)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl3,'linewidth',.5);
          ixi=or-.5;ix=ixi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-.5)*szx/length(ti)],[szy01 szy01+szy1]);set(h,'color',colorl3,'linewidth',.5);
          ixi=fi+.5;ix=fi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-or+.5)*szx/(fi-or+1)],[szy01 szy02+szy2]);set(h,'color',colorl3,'linewidth',.5);
          ixi=fi+.5;ix=ixi*DT;h=annotation('line',[szx0+(ixi-.5)*szx/length(ti) szx0+(ixi-.5)*szx/length(ti)],[szy01 szy01+szy1]);set(h,'color',colorl3,'linewidth',.5);          
      end
    
      
      h=plot([1:length(sigspk)]*DT,sigspk,'k+-');set(h,'linewidth',2);xlim(xl);set(gca,'fontsize',fs);ylim(yl);
  
    pause(.5);xt=get(gca,'xticklabel');xx=cell(size(xt,1),1);
    for ixi=1:length(xt),xx{ixi}=xt{ixi};end
    xx{round(length(xt)/2)}='time(ms)';
    set(gca,'xticklabel',xx);ylabel('{\Delta}F/F_0')
    h=line([xl(1) xl(1)],yl);set(h,'color','k');
    yt=get(gca,'ytick');
    for ixi=1:length(yt),
       h=line([xl(1) xl(1)+diff(xl)*.01],[yt(ixi) yt(ixi)]);set(h,'color','k'); 
    end
    
    axes('position',[szx0 szy02 szx szy2]);
    imagesc(catspk);set(gca,'fontsize',fs);set(gca,'xtick',[],'ytick',[]);
   
    %cm=jet(256);cm(1,:)=[0 0 0];colormap(cm)
    colormap(cm);
    % separadors blancs a la imatge dabaix
    for ixi=1:length(ti)
        ix=((ixi-1)*(size(catspk,1)+2))-0.5;
        h=line([ix,ix],[0,size(catspk,1)+1]);
        if(invertPatch==0),
            set(h,'color',colorp2,'linewidth',3);
        else
            set(h,'color',colora0,'linewidth',3);
        end
    end
    pos(4)=pos(4)*(szx*pos(3)/size(catspk,2))/(szy2*pos(4)/size(catspk,1));
    set(gcf,'color',colorf0,'position',pos);
    axes(HM);
          yt=get(gca,'ytick');
    for ixi=1:length(yt),
       h=line([xl(1) xl(1)+diff(xl)*.01],[yt(ixi) yt(ixi)]);set(h,'color','k'); 
    end
    set(gca,'ytick',yt);
    %remove text
    if(rt==1)
     set(gca,'ytick',[]);set(gca,'xtick',[]);ylabel('');
    end
      h=annotation('line',[szx0 szx0],[szy01 szy01+szy1]);set(h,'color','k');
      h=annotation('line',[szx0 szx0+szx],[szy01 szy01]);set(h,'color','k');



end