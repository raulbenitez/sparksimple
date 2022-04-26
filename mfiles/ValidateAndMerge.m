function [] = ValidateAndMerge()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escoger las carpetas que contienen los Results folder
% i darle a cancel para terminar.
%
% Puede que haga falta imponer un ResultsFolder tag si el ultimo
% SparkSimple se ha ejecutado con un ResultsFolder tag distinto 
% al que hay que validar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ch;
load('RF.mat');
RF='Results_170427';   % Descomentar esta linea si hace falta imponer un tag

fol=1;
N=0;
orpa='';
while(sum(fol)~=0)
    fol=uigetdir(orpa);
    if(sum(fol)~=0)
        orpa=fol;
        N=N+1;
        rutes{N}=fol;
    end
end
try,gg=urlread('http://leica.upc.es/counters/vm.php'); end
%rutes{1}='F:\SparkAnalysis\MAM13102301_CON_Pitx2_AD_s001';
RFN='MERGED';
scr=get(0,'ScreenSize');
taskbar=50;borders=7;
tamany=[scr(3)-2*borders scr(4)-1.5*taskbar-2*borders];%[1155 570];
posicio=[borders taskbar+borders ];%[round((scr(3)-tamany(1))/2) round((scr(4)-tamany(2))/2)];
% tamany=[1104         481];
% posicio=[1441         586 ];
options={'REJECT SPARK','ACCEPT SPARK'};
amp1=.4;
amp2=tamany(2)*0.43/tamany(1);
fp=10;fg=18;
altresu=.8;
if(1-3*fg/tamany(2)-3*fp/tamany(2)-altresu<.1),
    altresu=1-3*fg/tamany(2)-3*fp/tamany(2)-.1;
end
altsli=.04;

COL=[0 0 0];COL2=[0.4,0.2,0.2];
%COL=[.5 .1 .1];

gaga=0;
%%%%%%% test for good
failtest=0;
for ii=1:length(rutes) %
     
    S=dir([rutes{ii} '/' RF '*']);
    
    for jj=1:length(S) 
        if(S(jj).isdir==1)
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
%%%%%%%%%%%%% preparem GUI
W.main=figure(88);
set(W.main,'units','pixels','position',[posicio,tamany],...
    'menubar','none','name','ValidateAndMerge for SparkSimple.',...
    'numbertitle','off','resize','on','color',COL);
W.film=axes('position',[0 (1+altsli)/2 amp1 (1-altsli)/2],'xtick',[],'ytick',[],'color','k');
W.mask=axes('position',[0 0 amp1 (1-altsli)/2],'xtick',[],'ytick',[],'color','k');
W.slid=uicontrol('style','slider','units','normalized',...
    'position',[0 (1-altsli)/2 amp1 altsli],'min',0,'max',1,'SliderStep',[.01 .1],'value',1);
W.trac=axes('position',[amp1 0 amp2 1],'xtick',[],'ytick',[],'color','k');
W.tex1=uicontrol('style','text','unit','normalized',...
    'position',[amp1+amp2 1-3*fg/tamany(2) 1-amp1-amp2 3*fg/tamany(2)],'fontsize',fg,...
    'backgroundcolor',COL,'HorizontalAlignment','center',....
    'foregroundcolor','w','string','Initiating....');
W.tex2=uicontrol('style','text','unit','normalized',...
    'position',[amp1+amp2 1-3*fg/tamany(2)-3*fp/tamany(2) 1-amp1-amp2 3*fp/tamany(2)],'fontsize',fp,...
    'backgroundcolor',COL,'HorizontalAlignment','center',....
    'foregroundcolor','w','string','');
W.rads=uibuttongroup('units','normalized','pos',[amp1+amp2 altresu 1-amp1-amp2 1-3*fg/tamany(2)-3*fp/tamany(2)-altresu],...
    'backgroundcolor',COL2);
W.rad1=uicontrol(W.rads,'style','rad','unit','normalized',...
    'position',[.1 0 .4 1],'fontsize',fg,'min',0,'max',1,...
    'backgroundcolor',COL2,'string',options{1},'foregroundcolor','w',...
    'value',0,'HorizontalAlignment','center');
W.rad2=uicontrol(W.rads,'style','rad','unit','normalized',...
    'position',[.6 0 .4 1],'fontsize',fg,'min',0,'max',1,...
    'backgroundcolor',COL2,'string',options{2},'foregroundcolor','w',...
    'value',0,'HorizontalAlignment','center');
W.tex3=uicontrol('style','edit','unit','normalized',...
    'position',[amp1+amp2 0 1-amp1-amp2 altresu],'fontsize',fp,...
    'backgroundcolor',COL,'HorizontalAlignment','left','min',0,'max',2,....
    'foregroundcolor','w','string',['Use keys ''4'' and ''6'' to navigate in film and keys ''y'' and ''n'' to accept or reject.']);
%%%%%%%%%%%%%%%%%%
%%%%%% callbacks
set(W.main,'KeyPressFcn',{@keyPress1});set(W.slid,'KeyPressFcn',{@keyPress4});
set(W.tex1,'KeyPressFcn',{@keyPress6});set(W.tex2,'KeyPressFcn',{@keyPress7});
set(W.rad1,'KeyPressFcn',{@keyPress9});set(W.rad2,'KeyPressFcn',{@keyPress10});
set(W.tex3,'KeyPressFcn',{@keyPress11});

%%%%%
CMMMC=superjet(255);
roiR=2;%  microns

countSparks=0;
MASKS=0;clear MASKS;frameoff=0;
for ii=1:length(rutes) % recorrem carpetes
    ensenya(['folder ' num2str(ii) ': ' rutes{ii}]);
    
    S=dir([rutes{ii} '/' RF '*']);
    
    for jj=1:length(S)  % recorrem subcarpetes de resultats
        if(S(jj).isdir==1)
            
            % carreguem dades i inicialitzem GUI
            ensenya(['Reading files in ''' S(jj).name '''']);
            set(W.tex2,'string','Reading frames.');
            load([rutes{ii} '/' S(jj).name '/zMetaData.mat']);
            dst=round(roiR/DX);
            if(exist('MASKS','var')==0),MASKS=mask;else,try,MASKS=MASKS+mask;end;end
            mk=mask*.2;mk(:,:,2)=mask*.7;mk(:,:,3)=mk(:,:,1);
            axes(W.mask);imagesc(mk);axis equal;set(gca,'color','k');
            load([rutes{ii} '/' S(jj).name '/zData2.mat']);
            film=dir([rutes{ii} '/' S(jj).name '/detFilm/*.png']);
            [frame,map]=imread([rutes{ii} '/' S(jj).name '/detFilm/' film(1).name]);
            axes(W.film);imagesc(frame);axis equal;set(gca,'color','k');
            set(W.slid,'min',1,'max',length(film),'SliderStep',[1/length(film) 1/length(film)],'value',1);
            [a,b,c]=size(frame);
            peli=uint8(zeros(a,b,c,length(film)));peli(:,:,:,1)=frame;
            for kk=2:length(film)
                peli(:,:,:,kk)=imread([rutes{ii} '/' S(jj).name '/detFilm/' film(kk).name]);
            end
            set(W.tex2,'string','Reading raw images.');
            parrafada{1}='';cp=1;
            set(W.slid,'callback',{@slider,W,peli});
            load([rutes{ii} '/' S(jj).name '/zData0.mat']);
            vmax=max(max(max(volum)));
            fluoG=gray(255);fluoG(:,1)=0;fluoG(:,3)=0;
            nspk=length(spkF);
            % recorrem sparks
            set(W.tex3,'string','');
            ensenya('Manual validation');
            for kk=1:nspk
                axes(W.film);imagesc(peli(:,:,:,spkF(kk).pt));axis equal;set(gca,'color','k');
                set(W.slid,'value',spkF(kk).pt);
                set(W.tex1,'string',['Validate spark ' num2str(kk) '/' num2str(nspk) ' (frame' num2str(spkF(kk).pt) ')']);
                try
                if(spkF(kk).good==1),set(W.tex2,'string',['Accepted.']);set(W.rad2,'value',1);else set(W.tex2,'string',['Rejected. Reason: ' spkF(kk).fail]);set(W.rad1,'value',1);end
                catch
                    cprintf('string','Este error no debería ocurrir, contactar con alex y dadme el numero de experimento. Puede que sea un experimento procesado con una version antigua??');
                    error('spark structure has missing field ''good''');
                end
                [trace,map]=imread([rutes{ii} '/' S(jj).name '/spkFeat/Spk' num2str(kk) '.png']);
                map(1,:)=[0 0 0];map(2,:)=[1 1 1];
                axes(W.trac);imagesc(gray2rgb(trace,map));
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
                
                
                
                
                
                
                while(gaga==0) %%%%% esperem user input
                    pause();%%%%%%%%
                    on=round(get(W.slid,'value'));
                    if(ch=='4')&&(on>1),
                        set(W.slid,'value',on-1);
                        axes(W.film);imagesc(peli(:,:,:,on-1));axis equal;set(gca,'color','k');drawnow;
                        axes(W.mask);imagesc(volum(:,:,on-1),[0 vmax]);colormap(fluoG);axis equal;set(gca,'color','k');drawnow;
                    end
                    if(ch=='6')&&(on<size(peli,4))
                        set(W.slid,'value',on+1);
                        axes(W.film);imagesc(peli(:,:,:,on+1));axis equal;set(gca,'color','k');drawnow;
                        axes(W.mask);imagesc(volum(:,:,on+1),[0 vmax]);colormap(fluoG);axis equal;set(gca,'color','k');drawnow;
                    end
                    if(ch=='y')
                        ac=get(W.rad2,'value');
                        if(ac==1)
                            if(confirmed==1)
                                break;
                            else
                                confirmed=1;
                                set(W.tex1,'string',['Confirm accept?']);
                            end
                        else
                            set(W.tex1,'string',['Validate spark ' num2str(kk) '/' num2str(nspk) ' (frame' num2str(spkF(kk).pt) ')']);
                            set(W.rad2,'value',1) ;set(W.tex2,'string',['Change of mind']);confirmed=0;
                        end
                    end
                    if(ch=='n')
                        ac=get(W.rad1,'value');
                        if(ac==1)
                            if(confirmed==1)
                                break;
                            else
                                confirmed=1;
                                set(W.tex1,'string',['Negate again to confirm reject.']);
                            end
                        else
                            set(W.tex1,'string',['Validate spark ' num2str(kk) '/' num2str(nspk) ' (frame' num2str(spkF(kk).pt) ')']);
                            set(W.rad1,'value',1);set(W.tex2,'string',['Change of mind']);confirmed=0;
                        end
                    end
                end
                %textot=get(W.tex3,'string');
                if(ch=='y'),val=1;end
                if(ch=='n'),val=0;end
                parrafada=get(W.tex3,'string');cp=length(parrafada)+1;
                parrafada{cp}=['    Spark ' num2str(kk) ' : ' num2str(val)];
                set(W.tex3,'string',parrafada);
            end
            
            set(W.tex1,'string',['Folder finished. Confirm validation?']);
            gaga=0;confirmed=0;
            while(gaga==0)
                pause();
                if(ch=='y')
                    if(confirmed==1)
                        break;
                    else
                        confirmed=1;
                        set(W.tex1,'string',['Last chance to make changes in this folder. Move on to next folder?']);
                    end
                end
            end
            parrafada=get(W.tex3,'string');
            save([rutes{ii} '/' S(jj).name '/zValidation.mat'],'parrafada');
            set(W.tex1,'string','Please wait...');
            set(W.tex2,'string','Processing validation.');
            set(W.tex3,'string','');
            %%% afegim a la nova estructura de sparks
            for kk=1:length(parrafada)
                try
                    s1=strfind(parrafada{kk},'Spark');
                    s2=strfind(parrafada{kk},':');
                    s3=length(parrafada{kk});
                    spID=str2double(parrafada{kk}(s1+6:s2-1));
                    spGB=str2double(parrafada{kk}(s2+1:s3));
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
                    ensenya(['Could not tell if spark ' num2str(kk) ' was accepted or not. Assuming rejection.'])
                end
            end
            
        end
        frameoff=frameoff+temps;
    end
end

%%% tenim una superestructura SPARKS per a totes les carpetes i subcarpetes
ROIs=0;clear ROIs;
if(exist([rutes{1} '/' RFN])==0),mkdir([rutes{1} '/' RFN]);end
try,
save([rutes{1} '/' RFN '/lastValidation.mat'],'SPARKS','rutes','RFN');
end
ensenya('Validation saved');
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
im=zeros(size(ima));
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
    im=zeros(size(ima));
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
        end
        aux=round(mean(cooords,1));
        y=aux(1);x=aux(2);
        
        T=escriuFrase(num2str(ii),12);[sv,sh]=size(T);
        
        try im(y-dst+1:y-dst+sv,x-dst+1:x-dst+sh)=T;end
        try im(y-dst:y+dst,x-dst)=1;end
        try im(y-dst:y+dst,x+dst)=1;end
        try im(y-dst,x-dst:x+dst)=1;end
        try im(y+dst,x-dst:x+dst)=1;end
        
        coroi(ii,:)=[y x];
        
        %ROIs(ii).signal=tros;
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
        ROIs(ii).id=ii;
        ROIs(ii).Nspk=length(qui);
        ROIs(ii).signal=[1 2 3];
        
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
    
    try,im(MM==1)=1;ima(MM==1)=1;end
    ima(im==1)=1;
    IM=ima;IM(:,:,2)=ima;IM(:,:,3)=im;
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

[dades] = guardaDades2(SPARKS,ROIs,rutes{1},[rutes{1} '/' RFN]);


set(W.tex1,'string','Finished!');
set(W.tex2,'string','Merge successful.');


end

function [] = slider(varargin)

W=varargin{3};
peli=varargin{4};

on=get(W.slid,'value');
axes(W.film);imagesc(peli(:,:,:,round(on)));axis equal;set(gca,'color','k');drawnow;

end

function [] = keyPress1(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress2(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress3(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress4(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress5(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress6(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress7(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress8(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress9(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress10(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end
function [] = keyPress11(varargin)
global ch
ch=varargin{1,2}.Character;ch=lower(ch);
end