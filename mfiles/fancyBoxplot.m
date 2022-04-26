function [] = fancyBoxplot(D,X,C)

hold on;

if(nargin<2),N=1;end
N=unique(X);
if(nargin<3),C=zeros(length(N),3)+.4;end
dst=min(diff(N));

ms=30;% tamany markers
lw=3;% gruix ralles
sz=.7*dst;% llargada ralla gran
fa=.68;% factor entre ralla mes gran i mes peke
fcL=.4;% factor de color ralles (<1 per enfosquir)
res=68; % nombre de bins per determinar solapament
amp=.3; % factor d'ample maxim dels punts solapats respecte llargada ralla gran
qq=[.1 .5 .9];% quantils a mostrar (hauria de ser simetric)

ed=linspace(min(D),max(D),res);
for jj=1:length(N)
    


x=N(jj);
dist=D(X==x);


CP=C(jj,:);
CL=C(jj,:)*fcL;CL(CL>1)=1;

hi=histc(dist,ed);hh=zeros(size(hi));
sep=sz*amp/max(hi);
for ii=1:length(dist)
    on=find(ed<=dist(ii),1,'last');hh(on)=hh(on)+1;
    %hh(on)/hi(on)
    pos=x+(hh(on)*sep-sep*hi(on)*.5)-sep*.5;
    h=plot(pos,dist(ii),'.');   set(h,'markersize',ms,'color',CP);
    h=plot(pos,dist(ii),'ko');set(h,'markersize',ms*.3);% borde negre
end
qs=quantile(dist,qq);
aux=sz*(fa+(1-fa)*linspace(0,1,ceil(length(qs)/2)));aux=[aux(1:end-1) aux(end:-1:1)];
for ii=1:length(qs)
    h=line([x-aux(ii)/2 x+aux(ii)/2],[qs(ii) qs(ii)]);
    set(h,'linewidth',lw,'color',CL);
end

end

xlim([min(X)-dst/2 max(X)+dst/2]);

end