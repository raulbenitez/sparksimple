function [I] = ROIplot(imsize,x,y,titol,xlabl,ylabl,xlims,ylims,colos,gr,coo)
%Generates image with a plot.
% ImOfPlot(imsize,y) will plot vector y into an image where
%    imsize is a 2 or 3 component vector with the desired image
%    size in pixels: [ysize,xsize]
%
% ImOfPlot(imsize,x,y) will plot vector y versus vector x
%
% ImOfPlot(imsize,x,y,title,xlab,ylab) will include a title and labels on axes.
%
% ImOfPlot(imsize,x,y,title,xlab,ylab,xlims,ylims) can be used to specify
%    the limits in each axis. If not specified the limits are
%    the max and min value respectively.

gris=.6;
if(size(y,2)>size(y,1)),y=y';end

marge=ceil(min(imsize)*.16);
I=ones(imsize(1)-marge,imsize(2)-marge,3);
if(strcmp(titol,'')~=1),I=ones(imsize(1)-2*marge,imsize(2)-marge,3);end
[a,b,c]=size(I);
FS=round(marge*0.48);
% baix=min(y);dalt=max(y);
% esqu=min(x);dret=max(x);

%calcula regio on pintar
x0=[1:b];
xR=linspace(xlims(1),xlims(2),b);
if(xR(1)<x(1)),
    ori=find(xR>=x(1),1,'first');xR=xR(ori:end);x0=x0(ori:end);
else
    ori=find(x>=xR(1),1,'first'); x=x(ori:end);y=y(ori:end,:);
end
if(xR(end)<x(end)),
    fin=find(x<=xR(end),1,'last');x=x(1:fin);y=y(1:fin,:);
else
    fin=find(xR<=x(end),1,'last');xR=xR(1:fin);x0=x0(1:fin);
end
if(length(x0)>length(x)),
    x00=zeros(size(x));
    for ii=1:length(x)
        [aux,id]=min(abs(xR-x(ii)));
        x00(ii)=id;
    end
    x0=x0(x00);
    xR=xR(x00);
else
    x00=zeros(size(x0));
    for ii=1:length(x0)
        [aux,id]=min(abs(x-xR(ii)));
        x00(ii)=id;
    end
    x=x(x00);
    y=y(x00,:);
    aux=zeros(size(coo,1),1);
    for ii=1:size(coo,1)
        aux(ii)=find(x00>=coo(ii,1),1,'first');
    end
    coo(:,1)=aux;
end

% calcula ticks i enxufa grid
[ytiks] = niceTicks(ylims);
yimtiks=round(1+a-a*(ytiks-ylims(1))/(ylims(2)-ylims(1)));
ytiks(yimtiks<1)=[];yimtiks(yimtiks<1)=[];
ytiks(yimtiks>a)=[];yimtiks(yimtiks>a)=[];
[xtiks] = niceTicks(xlims);
try,[xtiks] = niceTicks(xlims,[1+length(xtiks) 2+length(xtiks)*2]);end
ximtiks=round(b*(xtiks-xlims(1))/(xlims(2)-xlims(1)));
xtiks(ximtiks<1)=[];ximtiks(ximtiks<1)=[];
xtiks(ximtiks>b)=[];ximtiks(ximtiks>b)=[];
for jj=1:length(yimtiks)
    I(yimtiks(jj),1:end,:)=gris;
end
for jj=1:length(ximtiks)
    I(:,ximtiks(jj),:)=gris;
end

% pinta grafics %%%%%%%%%%%%%%%%
for jj=1%:size(y,2)
    yy=y(:,jj);
    yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
    %yy(yy<=0)=1;yy(yy>=a)=a;
    for kk=2:length(x0)
        [allx1,ally1] = puntsEnmig([x0(kk-1) x0(kk)],[yy(kk-1) yy(kk)]);
        for kk2=1:length(allx1)
            lay=1+a-ally1(kk2);
            lax=allx1(kk2);
            if(((lay>0)&&(lay<=a))&&((lax>0)&&(lax<=b)))
                for cs=1:3
                    try
                        I(lay:lay+gr(jj),lax:lax+gr(jj),cs)=colos(jj,cs);
                    catch
                        gg=0;
                    end
                end
            end
        end
    end
end
I=I(1:a,1:b,:);

I=fletxes(I,x0,y,ylims,coo,colos,0);
I=fletxes(I,x0,y,ylims,coo,colos,1);
I=I(1:a,1:b,:);
% enxufa Y ticks
for jj=1:length(yimtiks)
    I(yimtiks(jj),1:round(marge*.5),:)=0;
end

% enxufa X ticks i grid
for jj=1:length(ximtiks)
    I(a-round(marge*.5):a,ximtiks(jj),:)=0;
end

% marc (eixos)
I(1,:,:)=0;I(end,:,:)=0;
I(:,1,:)=0;I(:,end,:)=0;

% afegeix marge
I=cat(2,ones(a,marge,3),I);
I=cat(1,I,ones(marge,b+marge,3));

% enxufa text per Y tiks
for jj=1:length(yimtiks)
    Z=text2im(num2str(ytiks(jj)));
    Z=imresize(Z,FS/size(Z,1));Z(Z<0)=0;Z(Z>1)=1;Z=Z/max(max(Z));
    if(size(Z,2)<marge)
        oriy=round(yimtiks(jj)-size(Z,1)/2);
        finy=oriy+size(Z,1)-1;
        orix=marge-size(Z,2);
        finx=orix+size(Z,2)-1;
        if(oriy>0)
            for cs=1:3
                I(oriy:finy,orix:finx,cs)=1-Z;
            end
        end
    end
end

% enxufa text per X tiks
for jj=1:length(ximtiks)
    Z=text2im(num2str(xtiks(jj)));
    Z=imresize(Z,FS/size(Z,1));Z(Z<0)=0;Z(Z>1)=1;Z=Z/max(max(Z));
    orix=round(marge+ximtiks(jj)-size(Z,2)/2);
    finx=orix+size(Z,2)-1;
    oriy=a+1;
    finy=oriy+size(Z,1)-1;
    if(finx<b+marge)
        for cs=1:3
            I(oriy:finy,orix:finx,cs)=1-Z;
        end
    end
end

% enxufa titol si cal
if(strcmp(titol,'')~=1)
    I=cat(1,ones(marge,b+marge,3),I);
    Z=text2im(titol);Z=imresize(Z,(FS+5)/size(Z,1));Z(Z<0)=0;Z(Z>1)=1;Z=Z/max(max(Z));
    oriy=round((marge/2)-(size(Z,1)/2));
    orix=round((size(I,2)/2)-(size(Z,2)/2));
    for cs=1:3
        I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1,cs)=1-Z;
    end
end

% xlabel
if(strcmp(xlabl,'')~=1)
    Z=text2im(xlabl);Z=imresize(Z,(2+FS)/size(Z,1));Z(Z<0)=0;Z(Z>1)=1;Z=Z/max(max(Z));
    oriy=round(size(I,1)-(size(Z,1)));
    orix=round((size(I,2)/2)-(size(Z,2)/2));
    for cs=1:3
        I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1,cs)=1-Z;
    end
end

% ylabel
if(strcmp(ylabl,'')~=1)
    Z=text2im(ylabl);Z=imresize(Z,FS/size(Z,1));Z(Z<0)=0;Z(Z>1)=1;Z=Z/max(max(Z));
    Z=Z';Z=Z(end:-1:1,:);
    orix=2;
    oriy=round((size(I,1)/2)-(size(Z,2)/1));
    for cs=1:3
        I(oriy:oriy+size(Z,1)-1,orix:orix+size(Z,2)-1,cs)=1-Z;
    end
end

% espai al final
I=cat(2,I,ones(size(I,1),1,3));

% tests de tamany
if(size(I,1)<imsize(1))
    I=cat(1,I,ones(imsize(1)-size(I,1),size(I,2),3));
end
if(size(I,2)<imsize(2))
    I=cat(2,I,ones(size(I,1),imsize(2)-size(I,2),3));
end
if(size(I,1)>imsize(1))
    I=I(1:imsize(1),:,:);
end
if(size(I,2)>imsize(2))
    I=I(:,1:imsize(2),:);
end

%test de rang
I(I<0)=0;I(I>1)=1;

end


function [I]=fletxes(I,x0,y,ylims,coo,colos,gutbad)
qq=find(coo(:,2)==gutbad);
coo=coo(qq,:);
%colos=colos(qq,:);
[a,b,c]=size(I);
% marca spks
BK0=zeros(a,b);BK=BK0;
pref='spk';numd=ceil(log10(max(coo(:,3))));gruc=0;grups=cell(0);
% comprova solapament
for jj=1:size(coo,1)
    yy=y(:,1);
    yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
    tco=[1+a-yy(coo(jj,1)) x0(coo(jj,1))]+[3 3];
    cul=round(tco+[-1 1]*min([round(tco(1))*.5,round(a/5)]));
    if(cul(2)+40>b),
        cul=round(tco+[-1 -1]*min([round(tco(1))*.5,round(a/5)]));
    end
    BK=textIm(cul(2),cul(1),[pref niceDigits(jj,numd)],zeros(a,b));
    BK(BK>0)=jj;
    q=unique(BK0(BK==jj));
    if((length(q)==1)&&(q==0)||(isempty(q)))
        gruc=gruc+1;
        grups{gruc}=[jj];
        BK0(BK==jj)=jj;
    else
        for ii=1:length(grups)
            if(~isempty(find(grups{ii}==q(end), 1))),break;end
        end
        grups{ii}=[grups{ii} jj];
    end
end
% exriu text i pinta fletxa
for jij=1:length(grups)
    qui=grups{jij};jj=qui(1);
    yy=y(:,1);
    yy=round(a*(yy-ylims(1))/(ylims(2)-ylims(1)));
    tco=[1+a-yy(coo(jj,1)) x0(coo(jj,1))]+[3 3];
    cul=round(tco+[-1 1]*min([round(tco(1))*.5,round(a/5)]));
    orient='left';
    if(cul(2)+40>b),cul=round(tco+[-1 -1]*min([round(tco(1))*.5,round(a/5)]));orient='right';end
    I=pintaFletxa(I,cul,tco,colos(1+qq(jj),:),1);
    fra='';
    for kk=1:length(qui)
        fra=[fra niceDigits(coo(qui(kk),3),numd) ','];
    end
    fra=fra(1:end-1);%;
    I=textIm(cul(2),cul(1),[pref fra],I,'textcolor',colos(1+qq(jj),:),'horizo',orient);
end

end