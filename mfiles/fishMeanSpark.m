function [sigspk,sumspk,catspk,number,DT,DX] = fishMeanSpark(Rfol)

% Rfol='F:\SparkAnalysis\071804_confo04_s015\Results';
% Rfol='H:\ValidatedSparks\13042503_CGS+MRS+ISO_s005\Results';
load([Rfol '/zMetaData.mat']);
load([Rfol '/zData0.mat']);
load([Rfol '/zData2.mat']);

roi=4;% roi radius
dst=round(roi/DX);

frames=[5 12];% frames to show before and after
tail=40;% samples to add in default tail

sumspk=zeros(2*dst+1,2*dst+1,sum(frames)+1);
sigspk=zeros(length(spkF(1).timeinterval)+tail,1);
catspk=[];
cnt=0;
for sp=1:length(spkF)
    
 if(spkF(sp).good==1)
     
             y=spkF(sp).py;
         x=spkF(sp).px;
    ts=spkF(sp).timesignal;
    ti=spkF(sp).timeinterval;ti=[ti(1):ti(end)+tail];
    
    if(ti(end)<=length(ts))
    ts=ts(ti);
    sigspk=sigspk+ts;
    centre=find(ts==max(ts));centre=centre(1); 
 
     or=centre-frames(1);fi=centre+frames(2);
     minivol=zeros(2*dst+1,2*dst+1,fi-or+1);
     if(fi<=length(ti)),
             for jj=or:fi
         minivol(:,:,jj-or+1)=volum(y-dst:y+dst,x-dst:x+dst,ti(jj));
         sumspk(:,:,jj-or+1)=minivol(:,:,jj-or+1)/length(or:fi);
         
             end
        cnt=cnt+1;
        
        imcat=[];
        for jj=1:size(minivol,3)
           imcat=[imcat zeros(2*dst+1,2)  minivol(:,:,jj)]; 
             
        end
        imcat=imcat(:,3:end);
        imcat=(imcat-min(min(imcat)))/(max(max(imcat))-min(min(imcat)));
 if(isempty(catspk)),catspk=imcat;else catspk=catspk+imcat;end
     end
    end
 end
 
end
number=cnt;
%disp(['Added ' num2str(cnt) ' sparks.']);





end