


rang=1000:1300;L=round(length(rang/2));
fix=3149;

figure(10);clf;set(gcf,'position',[231   390   672   718]);
axes('position',[0 .5 1 .5]);imagesc(Im);axis equal;
colormap([zeros(256,1),linspace(0,1,256)',zeros(256,1)]);
hold on;
h=line([fix-L,fix-L],[rang(1) rang(end)]);set(h,'color','b','linewidth',2);
h=line([fix+L,fix+L],[rang(1) rang(end)]);set(h,'color','b','linewidth',2);
h=line([fix-L,fix+L],[rang(end) rang(end)]);set(h,'color','b','linewidth',2);
h=line([fix-L,fix+L],[rang(1) rang(1)]);set(h,'color','b','linewidth',2);
set(gca,'xtick',[],'ytick',[]); 
% figure(11);
% plot(Im(:,fix));axis tight;
% set(gcf,'position',[ 1325         896        1019         339]);
axes('position',[0 0 1 .5]);
imagesc(Im(rang,fix-L:fix+L));colormap([zeros(256,1),linspace(0,1,256)',zeros(256,1)]);
axis equal;hold on;
h=line([2 2],[1 length(rang)]);set(h,'color','b','linewidth',2);
h=line([2*L 2*L],[1 length(rang)]);set(h,'color','b','linewidth',2);
h=line([1 2*L+1],[2 2]);set(h,'color','b','linewidth',2);
h=line([1 2*L+1],[length(rang)-1 length(rang)-1]);set(h,'color','b','linewidth',2);
h=line([L+1 L+1],[1 length(rang)]);set(h,'color','r','linewidth',2);
set(gca,'xtick',[],'ytick',[]);


figure(13);clf;h=plot(Im(rang,fix));set(h,'color','r','linewidth',2);
axis tight;ylim([0 .5]);
set(gcf,'position',[ 1325         896        1019         339]);
L1=.41;
L2=.31;
L3=.21;
Xt=280;
h=line([0 rang(end)],[L1 L1]);set(h,'color','k','linestyle','--');set(gca,'xtick',[]);
h=text(Xt,L1+.01,'$$ l_{i-1} $$','Interpreter', 'Latex');
set(h,'horizontalalignment','left','verticalalignment','bottom','color','k','fontsize',12);
h=line([0 rang(end)],[L2 L2]);set(h,'color','k','linestyle','--');set(gca,'xtick',[]);
h=text(Xt,L2+.01,'$$ l_i $$','Interpreter', 'Latex');
set(h,'horizontalalignment','left','verticalalignment','bottom','color','k','fontsize',12);
h=line([0 rang(end)],[L3 L3]);set(h,'color','k','linestyle','--');set(gca,'xtick',[]);
h=text(Xt,L3+.01,'$$ l_{i+1} $$','Interpreter', 'Latex');
set(h,'horizontalalignment','left','verticalalignment','bottom','color','k','fontsize',12);



I=imread('f1.png');
figure(10);
clf;set(gcf,'position',[231   390   672   718]);
axes('position',[0 0 1 1]);
imagesc(I);
h=line([468 2],[176 373]);set(h,'color','b','linewidth',2);
h=line([566 671],[176 373]);set(h,'color','b','linewidth',2);