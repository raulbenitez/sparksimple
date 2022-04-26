function [] = ModifyDist2MembraneManual()
close all;clear all;
load('RF.mat');
RF='Results_new';   % name used for results folders in spark simple
prefix='xDataD2M';save('prefix.mat','prefix');

fol=uigetdir('');

fol(fol=='/')='\';
 figure(1234);
 h=title(['Loading data..']);set(h,'fontsize',18);
    drawnow;       
       load([fol '/zMetaData.mat']); 
       indsep=find(fol=='\',1,'last');
if(strcmp(fol(indsep+1:end),'MERGED')==1)
    S=dir([fol(1:indsep) RF '*']);
    if(isempty(S)),error('Cant find results folder next to merge folder. Either one was moved or SparkSimple saved under a different name.');end
    load([fol(1:indsep) S(1).name '/zData0.mat']);
    load([fol  '/zDataMERGE.mat']);spkF=SPARKS;
else
    try,
       load([fol  '/zData0.mat']);  
       load([fol  '/zData2.mat']);
    catch
       error('Cant find result files from SparkSimple.'); 
    end
end
   

       fot=sum(volum,3);
       fot=(fot-min(min(fot)))/(max(max(fot))-min(min(fot)));
       [a,b]=size(fot);
       foto=zeros(a,b,3);
       f0=zeros(a,b);
       gr=ceil(min([a,b])/100);
       foto(:,:,2)=fot;
       f1=foto;


       
       for jj=1:length(spkF)
       if(spkF(jj).good==1)
          
           x=spkF(jj).px;
           y=spkF(jj).py;
           
           f1=foto;f2=foto;
           %%% contorn mask
           MM=mask>0;MM=MM-imerode(MM,[1 1 1;1 1 1;1 1 1]);
           fb=f2(:,:,1);fb(MM==1)=1;f2(:,:,1)=fb;
           fb=f2(:,:,2);fb(MM==1)=0;f2(:,:,2)=fb;
           fb=f2(:,:,3);fb(MM==1)=0;f2(:,:,3)=fb;
           f2(y-1:y+1,x,1)=1;f2(y,x-1:x+1,1)=1;
           f2(y-1:y+1,x,2)=1;f2(y,x-1:x+1,2)=1;
           f2(y-1:y+1,x,3)=0;f2(y,x-1:x+1,3)=0;
           %%%
           f1(y-1:y+1,x,1)=1;f1(y,x-1:x+1,1)=1;
           f1(y-1:y+1,x,2)=1;f1(y,x-1:x+1,2)=1;
           f1(y-1:y+1,x,3)=0;f1(y,x-1:x+1,3)=0;
           
           d2m=spkF(jj).dist2memb;
           d2m0=d2m;
           f1=enxufaLlegenda3(f1,DX,[1 1 0]);
          imagesc(f2);axis image;hold on;

          %%% cercle amb old d2m
          %rad=round(d2m/DX);angs=[1:360]*3.1416/180;
          %lesx=x+rad*cos(angs);lesy=y+rad*sin(angs);plot(lesx,lesy,'y')
          %%%
           h=title(['Spark ' num2str(jj) ', d2m= ' num2str(0.01*round(100*d2m)) '\mum.']);set(h,'fontsize',18);
           h=xlabel('(click outside image to accept)');set(h,'fontsize',14);
            flag=1;
           while(flag==1)
               
          
           [xs,ys]=ginput(1);
           if(((xs<=b)&&(ys<=a))&&((xs>0)&&(ys>0)))
           clf;
           d2m=norm([xs-x,ys-y])*DX;
                      imagesc(f1);axis image;hold on;
           h=title(['Spark ' num2str(jj) ', d2m= ' num2str(0.01*round(100*d2m)) '\mum.']);set(h,'fontsize',18);
           h=xlabel('(click outside image to accept)');set(h,'fontsize',14);
 
           h=line([xs x],[ys y]);set(h,'color','w');plot(xs,ys,'w+');
           else
               break;
           end
           end
           disp(['Spark ' num2str(jj) ' d2m : ' num2str(0.01*round(100*d2m0)) ' >> ' num2str(0.01*round(100*d2m))]);
           spkF(jj).dist2memb=d2m;
           
       end
       end
       
       if(strcmp(fol(indsep+1:end),'MERGED')==1)
       SPARKS=spkF;
       save([fol  '/zDataMERGE.mat'],'SPARKS','ROIs');
       [dades] = guardaDades(SPARKS,ROIs,fol(1:indsep-1),fol);
       else
       save([fol '/zData2.mat'],'spkF','ROIs');
       fol(fol=='\')='/';
       [dades] = guardaDades(spkF,ROIs,fol(1:indsep-1),fol);
       end
  disp('Folder finished!')
  close all;
 
end