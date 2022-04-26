function [dades] = guardaDades(spkF,ROIs,fol2,RF)

load('prefix.mat');

fol2(fol2=='\')='/';
titol=[fol2(find(fol2=='/',1,'last')+1:end)];


load([RF '/zMetaData.mat']);
%save([Rfol '/zMetaData.mat'],'mask','temps','DX','DT',);
meta{1}=[',,,,,,,,,,,,,,,,,,Pixel size(um)'];
meta{2}=['Frame rate(ms)'];
meta{3}=['Total duration(s)'];
meta{4}=['Cell area(um&sup2)'];
meta{5}=['Cell Fo'];
meta{6}=['Spark candidates'];

meta1{1}=[',,,,,,,,,,,,,,,,,,' num2str(DX)];
meta1{2}=[ num2str(DT)];
meta1{3}=[ num2str(temps*DT/1000) ];
meta1{4}=[ num2str(numel(find(mask==1))*DX*DX) ];
%meta1{5}=[ num2str(fo) ];
gg=reshape(fo,1,numel(fo));
gg(mask==0)=[];
meta1{5}=[niceNums(mean(gg),1) ' &#177 ' niceNums(std(gg),1)];
meta1{6}=[ num2str(length(spkF))];



ti=zeros(length(spkF),1);
for ii=1:length(spkF)
    ti(ii)=spkF(ii).pt;
end
[aux,ord]=sort(ti);
dades=[];
dadesOut=[];
%coords=[];
c=0;
countouts=0;
reasons=cell(0);
for izi=1:length(spkF)
    ii=ord(izi);
    if(spkF(ii).good==1)
        try
            aux=[spkF(ii).experimentFolderID spkF(ii).experimentSubFolderID];
            aux2=spkF(ii).pT;
        catch
            aux=[1,1];aux2=spkF(ii).pt;
        end
        
        %%%%  ID  X Y t d2m d2s AMP amp bl ror t2p FDHM FWHM tau roi
        dades=[dades;ii spkF(ii).px spkF(ii).py spkF(ii).pt aux2 spkF(ii).dist2memb spkF(ii).dist2spk spkF(ii).AMP spkF(ii).amp spkF(ii).BL spkF(ii).ror spkF(ii).mass spkF(ii).t2p spkF(ii).FDHM spkF(ii).FWHM spkF(ii).tau aux spkF(ii).roi];
        
    else
        countouts=countouts+1;
        try
            aux2=spkF(ii).pT;
        catch
            aux2=spkF(ii).pt;
        end
        reasons{countouts}=spkF(ii).fail;
        dadesOut=[dadesOut;ii spkF(ii).px spkF(ii).py aux2 spkF(ii).pt spkF(ii).dist2memb spkF(ii).AMP spkF(ii).amp spkF(ii).BL spkF(ii).ror spkF(ii).mass spkF(ii).t2p spkF(ii).FDHM spkF(ii).FWHM spkF(ii).tau];
    end
end

means=mean(dades);
stds=std(dades);
extraline='mean values & stds,';
means(1:4)=0;stds(1:4)=0;

dadesR=[];headr2=[];time=[];
if(~isempty(ROIs))
    for ii=1:length(ROIs)
        
        
        %%%%  ID  X Y t d2m d2s AMP amp bl ror t2p FDHM FWHM tau Nspk
        dadesR=[dadesR;ii ROIs(ii).px ROIs(ii).py 0 ROIs(ii).dist2memb ROIs(ii).dist2closest ROIs(ii).AMP ROIs(ii).amp ROIs(ii).BL ROIs(ii).ror ROIs(ii).mass ROIs(ii).t2p ROIs(ii).FDHM ROIs(ii).FWHM ROIs(ii).tau ROIs(ii).Nspk];
        
        
    end
    
    time=zeros(length(ROIs(1).signal),1+length(ROIs));
    time(:,1)=DT*([0:length(ROIs(1).signal)-1]');
    headr2=['time(ms),'];
    for jj=1:length(ROIs)
        headr2=[headr2 'ROI' num2str(jj) ','];
        for ii=1:length(ROIs(1).signal)
            time(ii,jj+1)=ROIs(jj).signal(ii);
            
        end
    end
end
meta{7}=['Surviving sparks'];
meta{8}=['Spk Freq (spk/s)'];
meta{9}=['Spk Dens (spk/s/nm&sup2)'];

meta1{7}=[num2str(length(spkF)-countouts)];
meta1{8}=[num2str((length(spkF)-countouts)/(temps*DT/1000)) ];
meta1{9}=[num2str(round((length(spkF)-countouts)/(temps*DT/1000)/(numel(find(mask==1))*DX*DX/1e6))) ];

%ensenya([num2str(length(spkF)-countouts) ' survived out of ' num2str(length(spkF)) '.']);
header='ID,Xpix,Ypix,Tframe,dist2memb,dist2closest,AMP,amp,BL,RoR,mass,t2p,FDHM,FWHM,tau';
header2=[header ',#sparks'];
header='ID,Xpix,Ypix,Tframe,Tframe,dist2memb,dist2closest,AMP,amp,BL,RoR,mass,t2p,FDHM,FWHM,tau,fldr,subfldr,ROI';
ensenya('Exporting files.');
muntaCSV([RF '/' prefix '_' titol '.csv'],',',meta,meta1,extraline,[means;stds],header,dades,'Roi Data',header2,dadesR,' ','ROI time signal',headr2,time);
spkHtml(RF,meta,meta1,dades,dadesR,dadesOut,reasons);

meta{7}=['Rejected sparks'];
meta{8}=['Rejected spark frequency(spk/s)'];
meta{9}=['Rejected spark density(spk/s/um²)'];

meta1{7}=[num2str(countouts)];
meta1{8}=[num2str((countouts)/(temps*DT/1000)) ];
meta1{9}=[num2str((countouts)/(temps*DT/1000)/(numel(find(mask==1))*DX*DX)) ];

header='ID,Xpix,Ypix,Tframe,Tframe,dist2memb,AMP,amp,BL,RoR,mass,t2p,FDHM,FWHM,tau';
muntaCSV([RF '/zDataRejected_' titol '.csv'],',',meta,meta1,header,dadesOut);



end