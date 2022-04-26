function [] = MergeOnly()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NG=1; % Especificar número de grupos de merge a hacer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escoger las carpetas que contienen los Results folder
% y darle a cancel para terminar cada grupo.
%
% Puede que haga falta imponer un ResultsFolder tag si el último
% SparkSimple se ha ejecutado con un ResultsFolder tag distinto
% al que hay que validar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RF.mat');
RF='Results_new';   % Descomentar esta linea si hace falta imponer un tag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orpa='';try load('lastmergefol.mat');orpa=fol;end
for ii=1:NG
    N=0;
    fol=1;
    rutes='1';clear rutes;
    while(sum(fol)~=0)
        fol=uigetdir(orpa,['Select folder #' num2str(N+1) ' for group #' num2str(ii)]);
        if(sum(fol)~=0)
            orpa=fol;
            N=N+1;
            rutes{N}=fol;
            save('lastmergefol.mat','fol');
        end
    end
    %orpa(orpa=='/')='\';
    %orpa=orpa(1:find(orpa=='\',1,'last')-1);
    if(exist('rutes','var')==0),error(['Expecting folders for merge in group ' num2str(ii)]);end
    RUTES{ii}=rutes;
end
%rutes{1}='F:\SparkAnalysis\MAM13102301_CON_Pitx2_AD_s001';
RFN='MERGED';


for geg=1:NG
    rutes=RUTES{geg};
    try,
        gg=urlread('http://leica.upc.es/counters/mo.php');
    end
    %%%%%%% test for good
    failtest=0;aux1=1;aux2=[];
    for ii=1:length(rutes) %
        
        S=dir([rutes{ii} '/' RF '*']);
        
        for jj=1:length(S)
            if(S(jj).isdir==1)
                load([rutes{ii} '/' S(jj).name '/zMetaData.mat']);try aux2(aux1)=roiR;aux1=aux1+1;end
                load([rutes{ii} '/' S(jj).name '/zData2.mat']);
                try
                    testvar=spkF(1).good;
                catch
                    ruu=rutes{ii};
                    ruu(ruu=='\')='/';
                    cprintf('error',['FAIL: folder ''' ruu '/' S(jj).name '''\n']);
                    cprintf('error',['Please reprocess with latest version of SparkSimple.\n']);
                    failtest=1;
                end
            end
        end
    end
    if(isempty(S)),error(['Found no results folders (searching tag: ''' RF ''').']);end
    if(failtest==1),error('Validation interrupted.');end
    %%%%%%
    if(aux1==1), ensenya(['ROI radius not found (data processed with old version), assuming 2um.'],superjet(1,'a'));roiR=2;end
    if(numel(unique(aux2))>1), ensenya(['ROI radius changes within experiments, assuming 2um.'],superjet(1,'a'));roiR=2;end
    % roiR=2;%  microns per minivol
    CMMMC=superjet(255);
    countSparks=0;
    MASKS=0;clear MASKS;frameoff=0;
    VOLUM=[];
    for ii=1:length(rutes) % recorrem carpetes
        ensenya(['folder ' num2str(ii) ': ' rutes{ii}]);
        
        S=dir([rutes{ii} '/' RF '*']);
        
        for jj=1:length(S)  % recorrem subcarpetes de resultats
            if(S(jj).isdir==1)
                countgood=0;
                % carreguem dades i inicialitzem GUI
                ensenya(['Reading files in ''' S(jj).name '''']);
                load([rutes{ii} '/' S(jj).name '/zMetaData.mat']);
                dst=round(roiR/DX);
                if(exist('MASKS','var')==0),MASKS=mask;else,try,MASKS=MASKS+mask;end;end
                % mk=mask*.2;mk(:,:,2)=mask*.7;mk(:,:,3)=mk(:,:,1);
                load([rutes{ii} '/' S(jj).name '/zData2.mat']);
                film=dir([rutes{ii} '/' S(jj).name '/detFilm/*.png']);
                [frame,map]=imread([rutes{ii} '/' S(jj).name '/detFilm/' film(1).name]);
                [a,b,c]=size(frame);
                
                load([rutes{ii} '/' S(jj).name '/zData0.mat']);
                try VOLUM=cat(3,VOLUM,volum);catch error('Images have different sizes!');end
                fluoG=gray(255);fluoG(:,1)=0;fluoG(:,3)=0;
                nspk=length(spkF);
                % recorrem sparks
                for kk=1:nspk
                    
                    gaga=0;confirmed=0;
                    
                    
                    %%%%% crea minivol de nou
                    
                    y=spkF(kk).py;
                    x=spkF(kk).px;
                    rgy6=max([1,y-dst]):min([a,y+dst]);
                    rgx6=max([1,x-dst]):min([b,x+dst]);
                    minivol=zeros(numel(rgy6),numel(rgx6),length(spkF(kk).timeinterval));
                    for ll=1:length(spkF(kk).timeinterval)
                        minivol(:,:,ll)=volum(rgy6,rgx6,spkF(kk).timeinterval(ll));
                        %          sumspk=sumspk+minivol(:,:,jj);
                    end
                    imcat=[];
                    for ll=1:length(spkF(kk).timeinterval)
                        imcat=[imcat zeros(numel(rgy6),1)  minivol(:,:,ll)];
                    end
                    imcat=(imcat-min(min(imcat)))/(max(max(imcat))-min(min(imcat)));
                    imwrite(uint8(255*imcat),CMMMC,[rutes{ii} '/' S(jj).name '/extraFigs/spk' num2str(kk) '.png']);
                    
                    
                    
                    
                    
                    
                    %%% afegim a la nova estructura de sparks
                    
                    try
                        
                        spID=kk;
                        spGB=spkF(spID).good;
                        countgood=countgood+spGB;
                        countSparks=countSparks+1;
                        SPARKS(countSparks).px=spkF(spID).px;
                        SPARKS(countSparks).py=spkF(spID).py;
                        SPARKS(countSparks).pt=spkF(spID).pt;
                        SPARKS(countSparks).pT=spkF(spID).pt+frameoff;
                        SPARKS(countSparks).FWHM=spkF(spID).FWHM;
                        SPARKS(countSparks).timesignal=spkF(spID).timesignal;
                        SPARKS(countSparks).timeinterval=spkF(spID).timeinterval;
                        SPARKS(countSparks).r2=spkF(spID).r2;
                        SPARKS(countSparks).amp=spkF(spID).amp;
                        SPARKS(countSparks).tau=spkF(spID).tau;
                        SPARKS(countSparks).ror=spkF(spID).ror;
                        SPARKS(countSparks).t2p=spkF(spID).t2p;
                        SPARKS(countSparks).FDHM=spkF(spID).FDHM;
                        SPARKS(countSparks).AMP=spkF(spID).AMP;
                        SPARKS(countSparks).BL=spkF(spID).BL;
                        SPARKS(countSparks).dist2memb=spkF(spID).dist2memb;
                        SPARKS(countSparks).mass=spkF(spID).mass;
                        SPARKS(countSparks).id=countSparks;
                        SPARKS(countSparks).good=spGB;
                        SPARKS(countSparks).fail=spkF(spID).fail;
                        SPARKS(countSparks).dist2spk=[];
                        SPARKS(countSparks).roi=[];
                        SPARKS(countSparks).experimentFolder=rutes{ii};
                        SPARKS(countSparks).experimentFolderID=ii;
                        SPARKS(countSparks).experimentSubFolder=S(jj).name;
                        SPARKS(countSparks).experimentSubFolderID=jj;
                        SPARKS(countSparks).experimentSubFolderSparkID=spkF(spID).id;
                    catch me
                        ensenya(['Could not tell if spark ' num2str(kk) ' was accepted or not. Assuming rejection.'],superjet(1,'h'))
                    end
                end
                frameoff=frameoff+temps;
            end
            
        end
    end
    
    %%% tenim una superestructura SPARKS per a totes les carpetes i subcarpetes
    ROIs=0;clear ROIs;
    if(exist([rutes{1} '/' RFN])==0),mkdir([rutes{1} '/' RFN]);end
    try,
        save([rutes{1} '/' RFN '/lastValidation.mat'],'SPARKS','rutes','RFN');
    end
    %ensenya('Validation saved');
    ensenya('Copying files');
    if(exist([rutes{1} '/' RFN '/detFilm'])==0),mkdir([rutes{1} '/' RFN '/detFilm']);end
    if(exist([rutes{1} '/' RFN '/extraFigs'])==0),mkdir([rutes{1} '/' RFN '/extraFigs']);end
    if(exist([rutes{1} '/' RFN '/spkFeat'])==0),mkdir([rutes{1} '/' RFN '/spkFeat']);end
    contaFrames=0;
    for ii=1:length(rutes) % recorrem carpetes per copiar frames peli
        S=dir([rutes{ii} '/' RF '*']);
        for jj=1:length(S)  % recorrem subcarpetes de resultats
            if(S(jj).isdir==1)
                SS=dir([rutes{ii} '/' S(jj).name '/detFilm/*.png']);
                for kk=1:length(SS)
                    contaFrames=contaFrames+1;
                    copyfile([rutes{ii} '/' S(jj).name '/detFilm/' SS(kk).name],[rutes{1} '/' RFN '/detFilm/F' num2str(100000+contaFrames) '.png']);
                end
            end
        end
    end
    cooords=[];spkind=[];
    ima=sum(volum,3);ima=(ima-min(min(ima)))/(max(max(ima))-min(min(ima)));
    ima0=ima;im=zeros(size(ima));
    ensenya('Merging spark structure.');
    for ii=1:length(SPARKS)  % recorrem sparks
        
        %  copia spark feature
        try
            ruta0=[SPARKS(ii).experimentFolder '/' SPARKS(ii).experimentSubFolder '/spkFeat/Spk' num2str(SPARKS(ii).experimentSubFolderSparkID) '.png'];
            ruta1=[rutes{1} '/' RFN '/spkFeat/Spk' num2str(ii) '.png'];
            copyfile(ruta0,ruta1);
        catch
            ensenya(['Missing image (Spk' num2str(SPARKS(ii).experimentSubFolderSparkID) ' in folder ' SPARKS(ii).experimentFolder '/' SPARKS(ii).experimentSubFolder ')']);
        end
        %  copia spark frames
        try
            ruta0=[SPARKS(ii).experimentFolder '/' SPARKS(ii).experimentSubFolder '/extraFigs/spk' num2str(SPARKS(ii).experimentSubFolderSparkID) '.png'];
            ruta1=[rutes{1} '/' RFN '/extraFigs/spk' num2str(ii) '.png'];
            copyfile(ruta0,ruta1);
        catch
            ensenya(['Missing image (Spk' num2str(SPARKS(ii).experimentSubFolderSparkID) ' in folder ' SPARKS(ii).experimentFolder '/' SPARKS(ii).experimentSubFolder ')']);
        end
        
        if(SPARKS(ii).good==1)
            y=SPARKS(ii).py;
            x=SPARKS(ii).px;
            spkind=[spkind ii];
            cooords=[cooords;y x];
            %im(y-1:y+1,x)=1;im(y,x-1:x+1)=1;
            im(y,x)=SPARKS(ii).id;
            
        end
    end
    %im=im(1:a,1:b);
    %ima(im==1)=1;
    imcr=im;imcr(imcr>1)=1;
    imcr=conv2(im,[0 1 0;1 1 1;0 1 0],'same');imcr(imcr>1)=1;
    im=labeler(im);im=im(1:a,1:b);
    try,
        MM=MASKS>0;MM=MM-imerode(MM,[1 1 1;1 1 1;1 1 1]);imcr(MM==1)=1;ima(MM==1)=1;
    end
    IM=imcr;IM(:,:,2)=ima;IM(:,:,3)=imcr;
    imwrite(IM,[rutes{1} '/' RFN '/extraFigs/LocationSparks2.png']);
    try,im(MM==1)=1;imac=ima;imac(im==1)=1;end
    IM=im;IM(:,:,2)=imac;IM(:,:,3)=im;
    imwrite(IM,[rutes{1} '/' RFN '/extraFigs/LocationSparks.png']);
    
    ensenya('Creating ROI structure.');
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
            SPARKS(spkind(ii)).dist2spk=min(aux);
        end
        coroi=zeros(length(unique(IDX)),2);
        Nrois=unique(IDX);
        for ii=1:length(Nrois)
            qui=find(IDX==ii);
            %%%%%
            buni(ii)=length(qui);
        end
        [ordRois ord]=sort(buni,'descend');
        
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
            for jj=1:length(qui)
                SPARKS(spkind(qui(jj))).roi=ii;
                y=SPARKS(spkind(qui(jj))).py;
                x=SPARKS(spkind(qui(jj))).px;
                cooords=[cooords;y x];
                fw=[fw SPARKS(spkind(qui(jj))).FWHM];
                r2=[r2 SPARKS(spkind(qui(jj))).r2];
                amp=[amp SPARKS(spkind(qui(jj))).FWHM];
                AMP=[AMP SPARKS(spkind(qui(jj))).AMP];
                tau=[tau SPARKS(spkind(qui(jj))).tau];
                ror=[ror SPARKS(spkind(qui(jj))).ror];
                t2p=[t2p SPARKS(spkind(qui(jj))).t2p];
                fd=[fd SPARKS(spkind(qui(jj))).FDHM];
                bl=[bl SPARKS(spkind(qui(jj))).BL];
                d2m=[d2m SPARKS(spkind(qui(jj))).dist2memb];
                mass=[mass SPARKS(spkind(qui(jj))).mass];
            end
            aux=round(mean(cooords,1));
            y=aux(1);x=aux(2);
            
            %        T=escriuFrase(num2str(ii),12);T(15:end,:)=[];[sv,sh]=size(T);
            y2=y;x2=x;mgh=ceil(roiR/DX);mgv=mgh;            try im(y-dst+1:y-dst+sv,x-dst+1:x-dst+sh)=T;end
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
            
            tros=VOLUM(oy:fy,ox:fx,1:end);
            
            tros=permute(mean(mean(tros)),[3 1 2]);
            ROIs(ii).signal=tros;
            ROIs(ii).px=x;
            ROIs(ii).py=y;
            ROIs(ii).Nsparks=jj;
            fw(isnan(fw))=[];ROIs(ii).FWHM=mean(fw);
            r2(isnan(r2))=[];ROIs(ii).r2=mean(r2);
            amp(isnan(amp))=[];ROIs(ii).amp=mean(amp);
            tau(isnan(tau))=[];ROIs(ii).tau=mean(tau);
            ror(isnan(ror))=[];ROIs(ii).ror=mean(ror);
            t2p(isnan(t2p))=[];ROIs(ii).t2p=mean(t2p);
            fd(isnan(fd))=[];ROIs(ii).FDHM=mean(fd);
            AMP(isnan(AMP))=[];ROIs(ii).AMP=mean(AMP);
            bl(isnan(bl))=[];ROIs(ii).BL=mean(bl);
            d2m(isnan(d2m))=[];ROIs(ii).dist2memb=mean(d2m);
            mass(isnan(mass))=[];ROIs(ii).mass=mean(mass);
            ROIs(ii).id=ii;
            ROIs(ii).Nspk=length(qui);
            
        end
        
        pd=squareform(pdist(coroi));
        if(size(coroi,1)==1),
            ROIs(1).dist2closest=0;
        else
            for ii=1:length(ROIs)
                aux=pd(ii,:);aux(aux==0)=[];
                ROIs(ii).dist2closest=min(aux)*DX;
            end
        end
        im=im(1:a,1:b,:);
        try,
            MM=mask>0;
            MM=MM-imerode(MM,[1 1 1;1 1 1;1 1 1]);
            im(MM==1)=1;ima(MM==1)=1;
        end
        IM=enxufaLlegenda3(im,DX);
        
        imwrite(IM,[rutes{1} '/' RFN '/extraFigs/LocationRois.png']);
        
        
    else
        if(numel(spkind)==1),spkF(spkind).dist2spk=0;spkF(spkind).roi=1;end
        ROIs=[];
    end
    
    
    try,
        mask=MASKS;mask(mask>0)=1;
    end
    temps=contaFrames;
    save([rutes{1} '/' RFN '/zMetaData.mat'],'mask','temps','DX','DT','fo');
    save([rutes{1} '/' RFN '/zDataMERGE.mat'],'SPARKS','ROIs');
    
    ROIplots(ROIs,SPARKS,DT,DX,roiR,[rutes{1} '/' RFN]);
    
    %[dades] = guardaDades2(SPARKS,ROIs,rutes{1},[rutes{1} '/' RFN]);
    [dades] = guardaDades(SPARKS,ROIs,rutes{1},[rutes{1} '/' RFN]);
    
    
    ensenya('Merge successful.');
    disp(' ');
    
    %cprintf('Strings',['Merge successful.\n']);
end
end
