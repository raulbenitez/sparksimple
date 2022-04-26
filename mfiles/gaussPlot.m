function [] = gaussPlot(D,X,C)


hold on;




if(nargin<2),N=1;X=ones(length(D),1);end
 N=unique(X);
if(nargin<3),list='k5cFA.geyorz';inds=randperm(length(list));C=superjet(length(N),list(inds(1:length(N))));end
if((nargin<2)||(length(N)==1)),dst=1;else dst=min(diff(N));end

% ms=30;% tamany markers
% lw=3;% gruix ralles
ext=.2; % espai extra per no tallar les distribucions
 sz=.7*dst;% amplitud curves
% fa=.68;% factor entre ralla mes gran i mes peke
 fcL=.4;% factor de color ralles (<1 per enfosquir)
 res=1000; % nombre de punts per pintar curves
% amp=.3; % factor d'ample maxim dels punts solapats respecte llargada ralla gran
qq=[.01 .5 .09];% quantils a mostrar (hauria de ser simetric)

rg=[min(D),max(D)];
ed=linspace(rg(1)-diff(rg)*ext,rg(2)+diff(rg)*ext,res);
for jj=1:length(N)
    


x=N(jj);
dist=D(X==x);

 [eix,curv]=roundBox(dist,res,[ed(1),ed(end)],qq);
%h=line([N(jj) N(jj)],[eix(end) eix(1)]);set(h,'Color',C(jj,:)*fcL); 
h=patch([N(jj)-curv*sz*.5,N(jj),N(jj)],[eix eix(end) eix(1)],C(jj,:));set(h,'FaceAlpha',.5,'EdgeColor',C(jj,:)*fcL);
h=patch([N(jj)+curv*sz*.5,N(jj),N(jj)],[eix eix(end) eix(1)],C(jj,:));set(h,'FaceAlpha',.5,'EdgeColor',C(jj,:)*fcL);

end
ylim([ed(1),ed(end)]);
xlim([min(X)-dst/2 max(X)+dst/2]);

set(gca,'box','on');grid on;

end

function [x,ex] = roundBox(dist,punts,rg,qtt)
% Retorna una curva composada de dues meitats de gaussiana centrades a la
% mitja i amb sigma igual al percentil 25/75. Normalitzada a 1.
% dist es el conjunt de valors 
% punts son el nombre de samples, rg son els extrems de x i qtt es per 
% especificar percentils diferents.

if(nargin<4),qtt=[.25 .50 .75];end
if(nargin<2),punts=100;end

y = quantile(dist,qtt);

y=[min(dist) y max(dist)];
if(nargin<3)
res=(max(y)-min(y))/(punts-1);
x=min(y):res:max(y);
else
res=(rg(2)-rg(1))/(punts-1);    
x=rg(1):res:rg(2);
end


A=1;
sig=y(4)-y(3);
ep=(x-y(3))/sig;ep=ep.*ep;
e1=A*exp(-ep/2);
sig=y(3)-y(2);
ep=(x-y(3))/sig;ep=ep.*ep;
e2=A*exp(-ep/2);

f=find(x>=y(3),1,'first');
ex=[e2(1:f-1),e1(f:end)];

% figure(1);clf;
% bar([0:.05:1],histc(dist,[0:.05:1]));
% hold on;yl=ylim();cc=0;xlim([-.05,1.05]);
% cc=cc+1;plot([y(cc) y(cc)],yl,'r');
% cc=cc+1;plot([y(cc) y(cc)],yl,'r');
% cc=cc+1;plot([y(cc) y(cc)],yl,'r','linewidth',2);
% cc=cc+1;plot([y(cc) y(cc)],yl,'r');
% cc=cc+1;plot([y(cc) y(cc)],yl,'r');
% plot(x,ex,'g','linewidth',2);

end