function [] = clusterROIs(fol,Rfol,SS,sr,bf)
load([Rfol '/zData0.mat']);
load([Rfol '/zData2.mat']);
load([Rfol '/zMetaData.mat']);
if(exist('ROIs','var')),clear ROIs;end
[a,b,c]=size(volum);
if(isempty(spkF)==0)
    ensenya('ROI clustering.');
    ima=sum(volum,3);ima=(ima-min(min(ima)))/(max(max(ima))-min(min(ima)));
    ima0=ima;
    im=zeros(size(ima));
    hema=ima0;hema(:,:,2)=0;hema(:,:,3)=ima0;hema(:,:,1)=0;
    CMMMC=jet(255);
    dst=round(roiR/DX);
    cooords=[];spkind=[];
    % np=zeros(length(spkF),2);tops=[a,b];qui=1;
    
    espuntdata=0;
    if ~strcmp(SS(1).name(end-2:end),'ata')
        Im=imread([fol '/' SS(1).name]);
    else
                xmlfil=dir([fol '/*.xml']);
        try
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(fol,xmlfil(1)) ;
        catch err
            [DX DY DT TTIME NXPIX NYPIX] = physical_parametersC(fol,xmlfil(1)) ;
        end
        
        Im=importdata([fol '/' SS(1).name]);
    if(size(Im,1)>1)&&(size(Im,2)>1)
        if((size(Im,1)==NYPIX)&&(size(Im,1))==NXPIX)            
        else
         NYPIX=size(Im,1);NXPIX=size(Im,2);   
        end
    else
        try
           Im=reshape(Im,[NYPIX NXPIX]); 
        catch           
     Im=reshape(Im,[150 52]);NYPIX=size(Im,1);NXPIX=size(Im,2);             
        end
    end
    end
    
    
    if(size(Im,3)==1)
        Ch=1;
    else
        [v,Ch]=max([sum(sum(Im(:,:,1))),sum(sum(Im(:,:,2))),sum(sum(Im(:,:,3)))]);
    end
    
    
    for ii=1:length(spkF)
        
        
        
        if(spkF(ii).good==1)
            
            
            cox=[spkF(ii).py,spkF(ii).px];
            while(im(cox(1),cox(2))>0)
                aux=1+round(rand(1));
                cox(aux)=cox(aux)+(round(rand(1))*2-1);
            end
            im(cox(1),cox(2))=ii;
            %          stronz=[spkF(ii).py,spkF(ii).px];
            %          if(qui==1),qui=2;else qui=1;end
            %          df=[spkF(ii).py,spkF(ii).px]-[a/2,b/2];
            %          %[aux,qui]=max(abs(stronz-[a/2,b/2]));
            %          if(df(qui)<0),
            %          np(ii,qui)=1;np(ii,3-qui)=stronz(3-qui);
            %          else
            %          np(ii,qui)=tops(qui);np(ii,3-qui)=stronz(3-qui);
            %          end
        end
    end
    im2=im;im2(im2>0)=1;
    if(length(unique(im))>1)
        try
            im2=labeler2(im,mask);im2=im2(1:a,1:b);%im=im2;
        catch
            im2=labeler(im);im2=im2(1:a,1:b);
        end
    end
    
    for ii=1:length(spkF)
        y=spkF(ii).py;
        x=spkF(ii).px;
        
        
        
        if(spkF(ii).good==1)

            spkind=[spkind ii];
            cooords=[cooords;y x];
            
            try
                sumspk=zeros(2*dst+1,2*dst+1);
                minivol=zeros(2*dst+1,2*dst+1,length(spkF(ii).timeinterval));
                for jj=1:length(spkF(ii).timeinterval)
                    minivol(:,:,jj)=volum(y-dst:y+dst,x-dst:x+dst,spkF(ii).timeinterval(jj));
                    sumspk=sumspk+minivol(:,:,jj);
                end
                sumspk=sumspk/jj;
                imcat=[];
                for jj=1:length(spkF(ii).timeinterval)
                    imcat=[imcat zeros(2*dst+1,1)  minivol(:,:,jj)];
                end
                imcat=(imcat-min(min(imcat)))/(max(max(imcat))-min(min(imcat)));
                imwrite(uint8(255*imcat),CMMMC,[Rfol '/extraFigs/spk' num2str(ii) '.png']);
            end
            
        end
    end
    im=im(1:a,1:b);
    ima(im>0)=1;
    IM=conv2(im,ones(ceil(b/268),ceil(b/268)),'same');IM(:,:,2)=im;IM(:,:,3)=im;IM(:,:,Ch)=ima;
    imwrite(IM,[Rfol '/extraFigs/LocationSparks.png']);
    IM=im2;IM(:,:,2)=ima0;IM(:,:,3)=im2;
    imwrite(IM,[Rfol '/extraFigs/LocationSparks2.png']);
    % figure;imagesc(uint8(255*IM));hold on;axis equal;
    % h=plot(cooords(:,2),cooords(:,1),'ow');set(h,'markersize',10);
    if(size(cooords,1)>1)
        cooords=cooords*DX;
        pd=pdist(cooords);
        h=linkage(pd);%figure;dendrogram(h)
        IDX = cluster(h,'Cutoff',roiR,'Criterion','distance');
        im=cat(3,ima0,ima0);im=cat(3,im,zeros(size(ima0)));
        
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
            fw=[];
            r2=[];
            amp=[];
            AMP=[];
            tau=[];
            ror=[];
            t2p=[];
            fd=[];
            bl=[];
            d2m=[];
            mass=[];
            for jj=1:length(qui)% ROI features are mean of forming good sparks
                spkF(spkind(qui(jj))).roi=ii;
                y=spkF(spkind(qui(jj))).py;
                x=spkF(spkind(qui(jj))).px;
                cooords=[cooords;y x];
                fw=[fw spkF(spkind(qui(jj))).FWHM];
                r2=[r2 spkF(spkind(qui(jj))).r2];
                amp=[amp spkF(spkind(qui(jj))).amp];
                AMP=[AMP spkF(spkind(qui(jj))).AMP];
                tau=[tau spkF(spkind(qui(jj))).tau];
                ror=[ror spkF(spkind(qui(jj))).ror];
                t2p=[t2p spkF(spkind(qui(jj))).t2p];
                fd=[fd spkF(spkind(qui(jj))).FDHM];
                bl=[bl spkF(spkind(qui(jj))).BL];
                d2m=[d2m spkF(spkind(qui(jj))).dist2memb];
                try,mass=[mass spkF(spkind(qui(jj))).mass];end
            end
            aux=round(mean(cooords,1));
            y=aux(1);x=aux(2);
            
            %        T=escriuFrase(num2str(ii),12);T(15:end,:)=[];[sv,sh]=size(T);
            y2=y;x2=x;mgh=ceil(roiR/DX);mgv=mgh;
            try im(y2-mgv:y2+mgv,x2-mgh,:)=1;end
            try im(y2-mgv:y2+mgv,x2+mgh,:)=1;end
            try im(y2-mgv,x2-mgh:x2+mgh,:)=1;end
            try im(y2+mgv,x2-mgh:x2+mgh,:)=1;end
            
            im=textIm(x,y,num2str(ii),im,'blending','on',...
                'textcolor',[1 1 1],'horizontalalignment','center','verticalalignment','mid');
            %figure(33);imagesc(imm);
            
            
            coroi(ii,:)=[y x];
            
            oy=max(1,y-dst);
            fy=min(y+dst,a);
            ox=max(1,x-dst);
            fx=min(x+dst,b);
            
            tros=volum(oy:fy,ox:fx,1:end);
            
            tros=permute(mean(mean(tros)),[3 1 2]);
            ROIs(ii).signal=tros;
            ROIs(ii).px=x;
            ROIs(ii).py=y;
            ROIs(ii).Nsparks=jj;
            ROIs(ii).FWHM=mean(fw);
            ROIs(ii).r2=mean(r2);
            ROIs(ii).amp=mean(amp);
            ROIs(ii).tau=mean(tau);
            ROIs(ii).ror=mean(ror);
            ROIs(ii).t2p=mean(t2p);
            ROIs(ii).FDHM=mean(fd);
            ROIs(ii).AMP=mean(AMP);
            ROIs(ii).BL=mean(bl);
            ROIs(ii).dist2memb=mean(d2m);
            ROIs(ii).mass=mean(mass);
            ROIs(ii).id=ii;
            ROIs(ii).Nspk=length(qui);
            
        end
        
        pd=squareform(pdist(coroi));
        if(~isempty(pd)),
            for ii=1:length(ROIs)
                aux=pd(ii,:);aux(aux==0)=[];
                ROIs(ii).dist2closest=min(aux)*DX;
                
            end
        else
            for ii=1:length(ROIs)
                ROIs(ii).dist2closest=0;
            end
        end
        im=im(1:a,1:b,:);
        try,
            MM=mask>0;
            MM=MM-imerode(MM,[1 1 1;1 1 1;1 1 1]);
            im(MM==1)=1;ima(MM==1)=1;
        end
        IM=enxufaLlegenda3(im,DX);
        %         ima=ima0;
        %         ima(im==1)=1;
        %         IM=ima;IM(:,:,2)=ima;IM(:,:,3)=im;
        
        
        imwrite(IM,[Rfol '/extraFigs/LocationRois.png']);
    else
        if(numel(spkind)==1),spkF(spkind).dist2spk=0;spkF(spkind).roi=1;end
        ROIs=[];
    end
    
    save([Rfol '/zData2.mat'],'spkF','ROIs');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ROIplots(ROIs,spkF,DT,DX,roiR,Rfol)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a,b,temps]=size(volum);
    % incompatible amb versions anteriors a 2016b
    %volum=volum.*mask;
    for ii=1:temps% alternativa
       volum(:,:,ii)=volum(:,:,ii).*mask;
    end    
    % figure,ed=[1:max(max(max(volum)))];bar(ed,histc(volum(:),ed));
    ed=[1:max(max(max(volum)))];hi=histc(volum(:),ed);hi=cumsum(hi/sum(hi));
    if(numel(volum)/1e6>100),th=.999;else th=.99;end
    mx=find(hi>=th,1,'first');
    volum=(volum-min(min(min(volum))))/(mx-min(min(min(volum))));
    
    % volum=NormArray(volum,0,1-quantile(volum(:),.9));
    %volum=volum*256;
    %%%%%%%%%%
    coords=[];
    for ii=1:length(spkF)
        
        if(((sr==1)&&(spkF(ii).good==0))||(spkF(ii).good==1))
            rang=spkF(ii).pt;rang=rang+[-1 0 1];% fotogrames on es pintara a la peli (3)
            for jj=1:length(rang)
                coords=[coords;spkF(ii).py spkF(ii).px rang(jj) ii spkF(ii).good];
            end
        end
    end
    %coords=round(coords);
    
    %peli=zeros(size(volum,1),size(volum,2),size(volum,3));
    ensenya('Building Film. 000%');
    
    for ii=1:temps
        
        percentatge(ii,temps);
        quins=find(coords(:,3)==ii);
        %F=M;F(:,:,2)=volum(:,:,ii);
        F=zeros(a,b);
        F(:,:,2)=zeros(a,b);
        F(:,:,3)=zeros(a,b);
        F(:,:,Ch)=volum(:,:,ii);
        F=F*bf;F(F>1)=1;
        for jj=1:length(quins),
%             te=escriuFrase(num2str(coords(quins(jj),4)),9);
%             [t1,t2]=size(te);
            
            if(coords(quins(jj),5)==1)
             colo=[1,1,1];
            else
             colo=[1,0,0];   
            end
            F=textIm(coords(quins(jj),2),coords(quins(jj),1),[' ' num2str(coords(quins(jj),4))],F,'blending','on','colour',colo);

            try % creueta vermella
                F(coords(quins(jj),1),max(1,coords(quins(jj),2)-2):min(b,coords(quins(jj),2)+2),1)=1;
                F(max(1,coords(quins(jj),1)-2):min(a,coords(quins(jj),1)+2),coords(quins(jj),2),1)=1;
                % creueta blanca
                if(coords(quins(jj),5)==1)
                    F(coords(quins(jj),1),max(1,coords(quins(jj),2)-2):min(b,coords(quins(jj),2)+2),2)=1;
                    F(max(1,coords(quins(jj),1)-2):min(a,coords(quins(jj),1)+2),coords(quins(jj),2),2)=1;
                    F(coords(quins(jj),1),max(1,coords(quins(jj),2)-2):min(b,coords(quins(jj),2)+2),3)=1;
                    F(max(1,coords(quins(jj),1)-2):min(a,coords(quins(jj),1)+2),coords(quins(jj),2),3)=1;
                end
            end
            
        end
        
F=textIm(1,1,[' ' num2str(ii)],F,'blending','on','colour',[1 1 1],'verticalalignment','top');

    
        imwrite(uint8(F*255),[Rfol '/detFilm/F' num2str(ii+100000) '.png']);
    end
end
end


