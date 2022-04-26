function [Fo] = sparkDetect(fol2,Rfol,SS,dm,tm,roiR,fs,nc,filtime,ss,forcedbaseline)

optB=1;% 0 per F/Fo, 1 per (F-Fo)/Fo

[DX,DT]=getRes(fol2);
if(DX==0)||(DT==0)
    try
        xmlfil=dir([fol2 '/*.xml']);
        try
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(fol2,xmlfil(1)) ;
        catch err
            [DX DY DT TTIME NXPIX NYPIX] = physical_parametersC(fol2,xmlfil(1)) ;
        end
    catch err
        error('Cannot read pixel size and frame rate!');
    end
end
disp(['         Resolution: ' niceNums(DX,2) 'um x ' niceNums(DX,2) 'um x ' niceNums(DT,2) 'ms.']);
disp(['         ROI size: ',num2str(1+2*round(roiR/DX)),'pix x ',num2str(1+2*round(roiR/DX)),'pix.']);
%filtime=40;%ms de filter window
%w=ceil(filtime/DT);%
w=filtime;if(w==0),w=1;end

espuntdata=0;
if ~strcmp(SS(1).name(end-2:end),'ata')
    Im=imread([fol2 '/' SS(1).name]);
else
    Im=importdata([fol2 '/' SS(1).name]);
    if(size(Im,1)>1)&&(size(Im,2)>1)
        if((size(Im,1)==NYPIX)&&(size(Im,1))==NXPIX)            
        else
         NYPIX=size(Im,1);NXPIX=size(Im,2);   
         ensenya(['Warning: Image dimensions do not match ''NumberOfElements'' field in xml file. Physical parameters may be wrong.'],superjet(1,'a'));
        end
    else
        try
           Im=reshape(Im,[NYPIX NXPIX]); 
        catch
            try
     Im=reshape(Im,[150 52]);NYPIX=size(Im,1);NXPIX=size(Im,2); 
ensenya(['Warning: Image dimensions do not match ''NumberOfElements'' field in xml file.'],superjet(1,'a'));
ensenya(['         Using possible working dimensions of [150 42]. Physical parameters may be wrong.'],superjet(1,'a'));
            catch
error('Image dimensions do not match ''NumberOfElements'' field in xml file. Could not reconstruct image.') 
            end
        end
    end
    espuntdata=1;
end

if(size(Im,3)==1)
    Ch=1;
else
    [v,Ch]=max([sum(sum(Im(:,:,1))),sum(sum(Im(:,:,2))),sum(sum(Im(:,:,3)))]);
    if(Ch==3), ensenya(' warning: blue channel contains more info');end
end
Im=Im(:,:,Ch);
[a,b,c]=size(Im);
volum=zeros(a,b,length(SS));

ensenya('Building Volume. 000%');
temps=length(SS);
fs=round(fs/DX);
fs=fs+mod(fs,2)+1;
GF=gaussiana2d(fs);
for ii=1:temps
    percentatge(ii,temps);
    if(espuntdata==0)
        Im=imread([fol2 '/' SS(ii).name]);
    else
        Im=importdata([fol2 '/' SS(ii).name]);
        Im=reshape(Im,[NYPIX NXPIX]); 
    end
    Im=Im(:,:,Ch);
    if(fs>0)
        %Im=medfilt2(Im,[fs,fs]);
        Im=conv2(double(Im),GF,'same');
    end
    volum(:,:,ii)=Im;
    
end

%%%%%%%%%%%%%%%%%%%%% mask

%mk=sum(volum,3);mk=(mk-min(min(mk)))/(max(max(mk))-min(min(mk)));mk=im2bw(mk,graythresh(mk));
if(espuntdata==0)
    mk=cellMask4(sum(volum,3));
else
    mk=ones(size(volum(:,:,1)));
end
imwrite(mk,[Rfol '/extraFigs/mask.png']);

%%%%%%%%%%%%%%%%%% Negative intensity
volum=ss*volum;
if(ss<0),volum=volum-min(min(min(volum)));end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% baseline

if(sum(size(forcedbaseline))==2)
    ensenya('Estimating Baseline.');
    vk=volum;for ii=1:temps,vk(:,:,ii)=mk;end
    vk=volum.*vk;vk=reshape(vk,numel(vk),1);vk(vk==0)=[];
    %normF=[quantile(vk,.05) quantile(vk,.95)];
    %     ensenya('Estimating Baseline. 000%');
    %     [lesy,lesx]=find(mk>-1);% cutrada per no fer doble loop
    %     estimatedFo=ones(size(mk));
    %     aux=length(lesx);
    %     for ii=1:aux
    %         percentatge(ii,aux);
    %         estimatedFo(lesy(ii),lesx(ii))=quantile(volum(lesy(ii),lesx(ii),:),.05);
    %     end
    
    estimatedFo=quantile(vk,.05);
    if(optB==0)
        volum=(volum)./estimatedFo;
    else
        volum=(volum-estimatedFo)./estimatedFo;
    end
else
    if(optB==0)
        volum=(volum)./forcedbaseline;
    else
        % compatible amb R2016b onwards
        %volum=(volum-forcedbaseline)./forcedbaseline;
        % previous versions
        for ii=1:temps
            volum(:,:,ii)=(volum(:,:,ii)-forcedbaseline)./forcedbaseline;
        end
    end
end

%%%%%%%%%%%%%%%%%%

%fprintf(['\b\b\b\b']);fprintf(['\n']);

ensenya('Wavelet transform. 000%');
V2=zeros(size(volum));
if(espuntdata==1)
    scala=[ceil(12/DT):ceil(48/DT)];
else
    scala=[5:10];
end
for ii=1:a
    percentatge(ii,a);
    for jj=1:b
        % if(mask(ii,jj)==1)
        
        sig=volum(ii,jj,:);
        sig=permute(sig,[3 1 2]);
        
        sig=smooth(sig,w);
        sig=sig-mean(sig);
        C=cwt(sig,scala,'gaus2');
        %figure;imagesc(C);
        C=sum(C,1);
        V2(ii,jj,:)=permute(C,[1 3 2]);
        
        % end
    end
end





% %


%ensenya('Filtering.');
% testin'
%umbral=.55; % sobre 1 el q es considera rellevant de la cwt
V3=V2;clear V2;
t3=ceil(2*dm/DX);%dm um
w3=ceil(2*tm/DT);%tm ms
%sizlim=.7; % if long axis longer than this, region is covering two sparks;
%B=ones(t3,t3,1);B=B/sum(sum(sum(B)));
B=gaussiana2d(t3);for ii=2:w3,B(:,:,ii)=B(:,:,1);end
B=B/sum(sum(sum(B)));
if((size(B,1)>size(V3,1))||(size(B,2)>size(V3,2))),ensenya('Error: image size is smaller than distance in microns for merging (dm).','r');end
if~isempty(B)
    ensenya('Smoothing.');
    V3=conv3(V3,B);
end




%V3=V3/max(max(max(V3)));
[ax1 ax2]=find(mk==1);
aux6=reshape(V3(round(mean(ax1)),round(mean(ax2)),:),1,temps);

umbral=mad(aux6)*2.22;% herencia del noisefactor (=1.5)

opt=0;
try
    if(opt==0)
        V3(V3<umbral)=0;
        IX=RegMax3D(V3,nc);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if(opt==1)
        V3(V3<umbral)=0;
        ensenya('3D Localmax. 000%');
        L=bwlabeln(V3>umbral,26);
        area=cell2mat(struct2cell(regionprops(L,'Area')));
        [~,ord]=sort(area,'descend');
        prop=regionprops(L,'Centroid');
        minvol=3*ceil((1/DX)*(1/DX)*(50/DT));% volum minim del connected component
        tots=numel(find(area>minvol));
        if(tots>temps/10),tots=round(temps/10);end
        IX=zeros(size(V3));
        for ii=1:tots
            percentatge(ii,tots+1);
            aux=round([prop(ord(ii)).Centroid]);
            IX(aux(2),aux(1),aux(3))=max(V3(L==ord(ii)));
        end
        percentatge(1,1);
    end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if(opt==2)
        cands=100;
        %busquem umbral que deixa cands zones
        mx=max(max(max(V3)));
        rg=[mx*.5,mx];
        valsa=[-1 -1];
        while 1
            vals=[max(max(max(bwlabeln(V3>rg(1),26)))) max(max(max(bwlabeln(V3>rg(2),26))))];
            next=mean(rg);
            [~,q]=min(abs(vals-100));
            rg(3-q)=next;
            if(diff(rg)==0),break;end
            if((vals(1)==valsa(1))&&(vals(2)==valsa(2))),break;end
            valsa=vals;
        end
        L=bwlabeln(V3>rg(1),26);
        prop=regionprops(L,'Centroid');
        IX=zeros(size(V3));
        tots=length(prop);
        for ii=1:tots
            percentatge(ii,tots);
            aux=round([prop(ii).Centroid]);
            IX(aux(2),aux(1),aux(3))=max(V3(L==ii));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
catch
    ensenya('Cancelled: No points above threshold',superjet(1,'r'));
    IX=zeros(size(V33));
end
%proximityTest(IX,dm,tm,DX,DT);

% use non candidate pixels for better Fo estimation
%ensenya('Reestimating Fo.');
if(sum(size(forcedbaseline))==2)
    % recuperem volum original
    if(optB==0)
        volum=volum.*estimatedFo;
    else
        volum=volum.*estimatedFo;volum=volum+estimatedFo;
    end
    ensenya('Reestimating Baseline. 000%');
    % ly=78;lx=401;gg=volum(ly-ceil(roiR/DX/2):ly+ceil(roiR/DX/2),lx-ceil(roiR/DX/2):lx+ceil(roiR/DX/2),:);gg=mean(gg,1);gg=mean(gg,2);gg=reshape(gg,temps,1);
    % figure,plot(gg)
    vaux=volum;suptemp=ceil(68/DT);
    B=gaussiana2d(ceil(2/DX));for ii=2:suptemp,B(:,:,ii)=B(:,:,1);end
    B=B/sum(sum(sum(B)));
    vaux=conv3(vaux,B);
    
    vaux(V3>0)=-1;
    
    reestimatedFo=zeros(size(mk));
    for ii=1:a
        percentatge(ii,2*a);
        for jj=1:b
            vec=reshape(vaux(ii,jj,:),[temps,1]);
            vec(vec==-1)=[];
            % figure,bar([min(vec):max(vec)],histc(vec,[min(vec):max(vec)]))
            if(isempty(vec))
                reestimatedFo(ii,jj)=estimatedFo;
            else
                %reestimatedFo(ii,jj)=quantile(vec,.05);
                reestimatedFo(ii,jj)=min(vec)+std(vec)*.5;
            end
        end
    end    
    Fo=reestimatedFo;
    try
        ra=round(roiR/DX);ra=ra-mod(ra,2)+1;ker=ones(ra,ra);
        aux=mk-imerode(mk,ker);
        gg=find(aux>0);mkz=zeros(size(mk));
        Foo=Fo;tots=length(gg);
        for ii=1:tots
            percentatge(2*ii,2*tots);
            thismk=mkz;
            thismk(gg(ii))=1;thismk=conv2(thismk,ker,'same');
            Foo(gg(ii))=max(max(thismk.*Fo));
        end
    end
    Fo=Foo;
    if(optB==0)
        volum=(volum)./Fo;
    else
        % compatible amb R2016b onwards
        %volum=(volum-forcedbaseline)./forcedbaseline;
        % previous versions
        for ii=1:temps
            volum(:,:,ii)=(volum(:,:,ii)-Fo)./Fo;
        end
    end
else
    Fo=forcedbaseline;
end
FFo=Fo.*mk;
cmap=superjet(255);fa=.4;
ff=gray2rgb(uint8(255*NormArray(FFo)),cmap);

cb=cbar('sz',[size(FFo,1),52,3],'cm',cmap,'cv',linspace(min(min(FFo)),max(max(FFo)),255),'fa',fa,'bg',[ff(1,1,1) ff(1,1,2) ff(1,1,3)],'tc',[1 1 1],'bc',[1 1 1]);
aux=sum(cb,3);aux=aux==3;
ff=cat(2,ff,cb);
%gg=(mk-imerode(mk,ones(3,3)));ff(gg==1)=.5;ff((a*size(ff,2))+find(gg==1))=.5;ff(2*(a*size(ff,2))+find(gg==1))=.5;
imwrite(ff,[Rfol '/extraFigs/Fo.png']);

%aveure(volum,IX,IX);
if(exist([Rfol '/spkFeat'],'dir')~=7),mkdir([Rfol '/spkFeat']);end

% try
yesno=0;

save([Rfol '/zData0.mat'],'volum','yesno');
save([Rfol '/zData1.mat'],'IX');

temps=size(volum,3);fo=Fo;mask=mk;
save([Rfol '/zMetaData.mat'],'fo','mask','temps','DX','DT','roiR');

end