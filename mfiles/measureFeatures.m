function [spkF] = measureFeatures(volum,IP2,DX,DT,Rfol,roiR)
load('ys.mat');
normVal=0;maxVal=ys;
%normVal=1.0;%ceil(max(max(max(volum))));
delete([Rfol '/spkFeat/Spk*.*']);
margX=round(1*roiR/DX);
convs=round(1/DX);
if(mod(convs,2)==0),convs=convs-1;end
F=gaussiana2d(convs);
margT=round(160/DT);

sc=get(0,'screensize');
hori=ceil(.25*(sc(3)-60));
vert=ceil(sc(4)*3/5);

[a0,b0,c0]=size(IP2);

spk=numel(find(IP2>0));



IP=zeros(size(IP2));


vals=sort(unique(IP2),'descend');
%%%%%%%% opt old 5 color png
COLORS=[1 1 1;0 0 0;.8 .1 .1;.1 .8 .1;.1 .1 .8;.6 .6 .6];
%%%%%%%% opt new blending
COLORS2=superjet(6,'wkVji5');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contaspk=0;
ensenya('Spark Features. 00%');fprintf('\b');
for ii=1:length(vals)-1
    %ii=ii+1;
    thresh=vals(ii);
    
    cand=find(IP2==thresh);
    
    for jj=1:length(cand)
        contaspk=contaspk+1;
        
        numeret=round(100*contaspk/spk);
        if(numeret<10),numeret=['0' num2str(numeret)];else numeret=num2str(numeret);end
        fprintf(['\b\b\b\b ' numeret '%%']);
        
        [a1,b1,c1]=ind2sub([a0,b0,c0],cand(jj));
        
%         if(contaspk==22)
%             kkk=0;
%         end
        
        spkF(contaspk).px=b1;
        spkF(contaspk).py=a1;
        spkF(contaspk).pt=c1;
        
        oy=max(1,a1-margX);
        fy=min(a1+margX,a0);
        ox=max(1,b1-margX);
        fx=min(b1+margX,b0);
        ot=max(1,c1-margT);
        ft=min(c1+2*margT,c0);
        oy2=max(1,a1-round(1.5*margX));
        fy2=min(a1+round(1.5*margX),a0);
        ox2=max(1,b1-round(1.5*margX));
        fx2=min(b1+round(1.5*margX),b0);
        ot2=max(1,c1-round(margT/3));
        ft2=min(c1+round(2*margT/3),c0);
        ce=[a1-oy+1,b1-ox+1,c1-ot+1];
        
        tros0=volum(oy:fy,ox:fx,ot:ft);
        tr0=sum(volum(oy2:fy2,ox2:fx2,ot2:ft2),3);
        tr0=conv2(tr0,F,'same');
%         bg=mean([tr0(convs:end-convs,convs);tr0(convs:end-convs,end-convs);tr0(convs,convs:end-convs)';tr0(end-convs,convs:end-convs)']);
%         if(isnan(bg)),bg=min(min(tr0));end
        bg=min(min(tr0));
        tr0=tr0-bg;tr0(tr0<0)=0;
        tr0=(tr0)/(max(max(max(tr0))));
        [yy,xx]=find(tr0==max(max(tr0)));
        %tros1=IP2(max(1,a1-margX):min(a1+margX,a0),max(1,b1-margX):min(b1+margX,b0),max(1,c1-margT):min(c1+2*margT,c0));
        if(size(tr0,1)~=size(tr0,2)),extr=' candidate is too close to image limits.     ';else extr='.     ';end
        try
       % R=BlobRadiusC(tr0,[yy(1) xx(1) ce(3)]);
        R=DX*Blob4Diameters(tr0)/2;
        catch
            fprintf(['\nWarning: Could not measure spark radius (spk#' num2str(contaspk) ')' extr]);
            R=0;
        end
        spkF(contaspk).FWHM=2*R;
        
%         r=round(R/2);
%         oy=max(1,a1-r);
%         fy=min(a1+r,a0);
%         ox=max(1,b1-r);
%         fx=min(b1+r,b0);
        
        ts=permute(volum(oy:fy,ox:fx,1:c0),[3,1,2]);
        ts=mean(mean(ts,3),2);
        spkF(contaspk).timesignal=ts;
        spkF(contaspk).timeinterval=[ot:ft];
        cut=ts(ot:ft);
        normVal=max(normVal,max(cut));
        if(normVal>maxVal),normValNow=maxVal;else,normValNow=normVal;end
        %normVal=0;maxVal=ys;
        try
        [sortida r2 I] = sparkParametersFast2([vert hori],[ot:ft]*DT,cut,c1-ot,1,['Spk ' num2str(contaspk)],length(cut),normValNow,COLORS2);    
       % figure,imagesc(I);
        imwrite(I,[Rfol '/spkFeat/Spk' num2str(contaspk) '.png']);
        catch
            I=1;sortida=[NaN,NaN,NaN,NaN,NaN,NaN,NaN];r2=NaN;
            try
        [sortida r2 I] = sparkParametersFast2([1000 430]*2,[ot:ft]*DT,cut,c1-ot,1,['Spk ' num2str(contaspk)],length(cut),normValNow);
            end
        I(1,1)=5;
        imwrite(I+1,COLORS,[Rfol '/spkFeat/Spk' num2str(contaspk) '.png']);
        end
        % sortida=[amp,tau,ror,t2p,dhm,AMP];
        
        spkF(contaspk).r2=r2;
        spkF(contaspk).amp=sortida(1);
        spkF(contaspk).tau=sortida(2);
        spkF(contaspk).ror=sortida(3);
        spkF(contaspk).t2p=sortida(4);
        spkF(contaspk).FDHM=sortida(5);
        spkF(contaspk).AMP=sortida(6);
        spkF(contaspk).BL=sortida(6)-sortida(1);
        spkF(contaspk).mass=sortida(7);
        
    end
    
end

fprintf(['\b\b\b\b']);fprintf(['\n']);



if( exist('spkF','var')==0 ),spkF=1;end
    









end