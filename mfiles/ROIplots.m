function [] = ROIplots(ROIs,spkF,DT,DX,roiR,Rfol)
if(exist([Rfol '/roiFeat'],'dir')~=7),mkdir([Rfol '/roiFeat']);end
load('ys.mat');
sc=get(0,'screensize');
hori=sc(3)-60;
vert=ceil(hori*.1);
dist=roiR/DX;
for ii=1:length(ROIs)

    sig=ROIs(ii).signal;
    cenR=[ROIs(ii).py ROIs(ii).px];
    spksig=zeros(length(sig),1);coo=[];
    cc=1;
    cols=[0,0,0];ex=superjet(2,'ez');
    for jj=1:length(spkF)
    posS=[spkF(jj).py spkF(jj).px];
    aux=spkF(jj).roi;if(isempty(aux)),aux=0;dist=1.4*roiR/DX;else dist=roiR/DX;end
        if((norm(cenR-posS)<dist)||(aux==ii))
            try
           spksig(:,cc)=spkF(jj).timesignal;
            end
            try aux=spkF(jj).pT;catch aux=spkF(jj).pt;end
           coo(cc,1)=aux;
           if(spkF(jj).good==1)
               coo(cc,2)=1;
           cols(cc+1,1:3)=ex(1,:);
           else
               coo(cc,2)=0;
           cols(cc+1,1:3)=ex(2,:);
           end
           coo(cc,3)=jj;
           cc=cc+1;
        end
    end
    data=[sig,spksig];
    yl=[min(min(data)) max(max(data))];
    yl=[min(yl(1),0) max(ceil(yl(2)),1.5)];
    %spark arrows
    [I] = ROIplot([vert hori],[1:length(sig)]*DT,data,['ROI' num2str(ii)],'time [ms]','',[0,length(sig)*DT],yl,cols,[2,ones(1,size(spksig,2))],coo);
%   figure,imagesc(I),axis image
    imwrite(uint8(I*255),[Rfol '/roiFeat/ROI' num2str(ii) '.jpg']);
    
end

end