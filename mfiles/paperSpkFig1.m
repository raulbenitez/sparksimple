
%%  figureta paper sparks

clc;clear all;close all;
spk=2;
expfol='C:\Users\avallmit\Documents\OwnPapers\sparks\09031701-01-seq-s003_R-80_CON2\Results_new2\';


load([expfol 'zData2.mat']);
load([expfol 'zMetaData.mat']);


COLORS2=superjet(51,'wkrgb5');
COLORS2=COLORS2([1:10:51],:);
COLORS=[1 1 1;0 0 0;.8 .1 .1;.1 .8 .1;.1 .1 .8;.5 .5 .5];


ot=spkF(spk).timeinterval(1);
ft=spkF(spk).timeinterval(end)+23;
c1=spkF(spk).pt-ot;
ts=spkF(spk).timesignal;
cut=ts(ot:ft);
ft=ft-ot;ot=0;
normVal=1.7;



[sortida r2 I] = gggg([1000 430],[ot:ft]*DT,cut,c1,1,['Spark ' num2str(spk)],length(cut),normVal);

[sortida r2 I] = sparkParametersFast2([600 300],[ot:ft]*DT,cut,c1,1,['Spark ' num2str(spk)],length(cut),normVal,COLORS2);
imwrite(I,['Spk' num2str(spk) 'b.png']);
[sortida r2 I] = sparkParametersFast2([1000 430],[ot:ft]*DT,cut,c1,1,['Spk ' num2str(spk)],length(cut),normVal);
I(1,1)=5;
imwrite(I+1,COLORS,['Spk' num2str(spk) '.png']);


x=[ot:ft]*DT;
T=cut;
centre=c1;
% AMP
[amp ind]=max(T(centre(1)-3:centre(1)+8));% (F-Fo)/Fo
centre(1)=centre(1)-3+ind-1;
bline=quantile(T(1:centre(1)),.02);
sortida(6)=amp;
amp=amp-bline;
sortida(1)=amp;



aux1=centre(1);
uprise=T(1:aux1);
di=zeros(1,aux1-1);
for ii=1:aux1-1
    t1=uprise(1:ii);
    t2=uprise(ii:centre(1));
    [fresult2,S]=polyfit([ii:centre(1)]',t2,1);
    t3=[ii:centre(1)]*fresult2(1)+fresult2(2);
    %     figure;plot(uprise);hold on;
    %     plot([ii:centre(1)],t3,'r');
    %     plot([1:ii],mean(t1)*ones(1,ii),'r');
    ft=[mean(t1)*ones(1,ii) [ii:centre(1)]*fresult2(1)+fresult2(2)];
    ft(ii)=(ft(ii)+ft(ii+1))/2;
    ft(ii+1)=[];
    
    di(ii)=sum(abs(ft'-uprise));
end
[kk,aux]=min(di);

% %%%%%%%%%%%%%%
T1=T(aux:aux1);
x1=x(aux:aux1);
%figure;plot(x,T);hold on;plot(x1,T1,'r');

% fit lineal a la pujada
if(size(x1,1)<size(x1,2)),x1=x1';end
if(size(T1,1)<size(T1,2)),T1=T1';end
[fresult2,S]=polyfit(x1,T1,1);
T1=fresult2(1)*x1+fresult2(2);
fitted1=fresult2(1)*x+fresult2(2);
%figure;plot(x,fitted1);hold on;plot(x,fitted2);


% t2p
%aux=find(T(1:centre(1))<1,1,'last');
xb=(bline-fresult2(2))/fresult2(1);
xp=x(aux+round((centre(1)-aux)/2));
t2p=x(centre(1))-xb;
sortida(4)=t2p;


% RoR
ror=fresult2(1);
sortida(3)=ror;



% FDHM
aux=centre(1)-1+find(T(centre(1):end)<bline+amp/2,1,'first');
aux2=find(T(1:centre(1))<bline+amp/2,1,'last');
if((~isempty(aux))&&(~isempty(aux2))),dhm=x(aux)-x(aux2);x51=x(aux2);x52=x(aux);h50=T(aux); else dhm=NaN;end
sortida(5)=dhm;
% Tau
caiguda=T(centre(1):end);
TT=caiguda;
xx=x(centre(1):end);



% fit exponencial al decay
model1 = fittype('a*exp(-b*(x-e))+f');
options1 = fitoptions(model1);
options1.Display = 'off';
options1.Robust = 'bisquare';
%%%%%%%%%%%%%%%%%%%%%%      a         b        e        f
options1.startpoint = [max(abs(TT)) 0  min(xx) min(TT)];
options1.Lower = [0 0 -10*max(xx) min(TT)];
options1.Upper = [100*max(abs(TT)) 1 100*max(xx) max((TT))];
[fresult1,gof2,out2]  = fit(xx',TT,model1,options1);
rsquare1=gof2.rsquare;
r2=rsquare1;
%     figure;plot(x,T,'-b'); hold on;plot(fresult1,'r');
fitted2=fresult1.a*exp(-fresult1.b*(x-fresult1.e))+fresult1.f;
tau=1/fresult1.b;




figure;plot(x,T),hold on,plot(x,fitted1),plot(x,fitted2)
ylim([-.1,normVal]),

muntaCSV(['spk' num2str(spk) '.csv'],';','time(ms);spark;uprise;decay;',[x',T,fitted1',fitted2']);


%%

roiR=2;
load([expfol 'zData0.mat']);
ima=sum(volum,3);ima=(ima-min(min(ima)))/(max(max(ima))-min(min(ima)));
ima0=ima;
im=zeros(size(ima));
[a,b]=size(im);
hema=im;
CMMMC=jet(255);
dst=round(roiR/DX);
cooords=[];spkind=[];
% np=zeros(length(spkF),2);tops=[a,b];qui=1;
for ii=1:length(spkF)
    
    if(spkF(ii).good==1)
        cox=[spkF(ii).py,spkF(ii).px];
        while(im(cox(1),cox(2))>0)
            aux=1+round(rand(1));
            cox(aux)=cox(aux)+(round(rand(1))*2-1);
        end
        im(cox(1),cox(2))=ii;
        
    end
end
im2=im;im2(im2>0)=1;
if(length(unique(im))>1)
    im=labeler(im);
end


for ii=1:length(spkF)
    
    if(spkF(ii).good==1)
        %try
        y=spkF(ii).py;
        x=spkF(ii).px;
        spkind=[spkind ii];
        cooords=[cooords;y x];
        
    end
    
end
im=im(1:a,1:b);

ima(im==1)=1;
IM=im;IM(:,:,2)=ima;IM(:,:,3)=im;
imwrite(IM,['LocationSparks.png']);
IM=im2;IM(:,:,2)=ima0;IM(:,:,3)=im2;
imwrite(IM,['LocationSparks2.png']);

%%%%%%%
if(size(cooords,1)>1)
    cooords=cooords*DX;
    pd=pdist(cooords);
    h=linkage(pd);%figure;dendrogram(h)
    IDX = cluster(h,'Cutoff',roiR,'Criterion','distance');
    im=zeros(size(ima));
    %dst=dst*2;
    
    
    sq=squareform(pd);
    for ii=1:length(spkind)
        aux=sq(ii,:);aux(aux==0)=[];
        spkF(spkind(ii)).dist2spk=min(aux);
    end
    
    
    coroi=zeros(length(unique(IDX)),2);
    Nrois=unique(IDX);
    for ii=1:length(Nrois)
        qui=find(IDX==ii);
        %%%%%
        buni(ii)=length(qui);
    end
    [ordRois ord]=sort(buni,'descend');
    %%%%%
    for ii=1:length(Nrois)
        qui=find(IDX==ord(ii));
        cooords=[];
        for jj=1:length(qui)
            
            y=spkF(spkind(qui(jj))).py;
            x=spkF(spkind(qui(jj))).px;
            cooords=[cooords;y x];
            
        end
        aux=round(mean(cooords,1));
        y=aux(1);x=aux(2);
        T=escriuFrase(num2str(ii),12);T(15:end,:)=[];[sv,sh]=size(T);
        y2=y;x2=x;mgh=round(sh/2);mgv=round(sv/2);
        try im(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
        try im(y2-mgv-1:y2-mgv+sv,x2-mgh-1)=1;end
        try im(y2-mgv-1:y2-mgv+sv,x2-mgh+sh+1)=1;end
        try im(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
        try im(y2-mgv+sv+1,x2-mgh-1:x2-mgh+sh+1)=1;end
     %   im=textIm(x,y,num2str(ii),im,'horizontalalignment','center','verticalalignment','mid','blending','off','background','circle','textcolor',1,'bdcolor',1,'bgcolor',0);
        
        
        
    end
    
    try,
        MM=mask>0;
        MM=MM-imerode(MM,[1 1 1;1 1 1;1 1 1]);
        im(MM==1)=1;ima(MM==1)=1;
    end
    im=im(1:a,1:b);
    ima=ima0;
    ima(im==1)=1;
    IM=ima;IM(:,:,2)=ima;IM(:,:,3)=im;
    
    
    imwrite(IM,['LocationRois.png']);
end