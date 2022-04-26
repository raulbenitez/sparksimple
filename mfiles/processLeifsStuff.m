function [] = processLeifsStuff()

ii=0;
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0005_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0006_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0007_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0009_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0010_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0011_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0017_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0019_2808_mask-wrong';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0020_2808_mask-wrong';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0021_2808_mask-wrong';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0023_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0024_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0025_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0026_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0028_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0030_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0031_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0049_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0052_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif/HAM-0053_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0054_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0059_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0079_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0084_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0115_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0116_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0120_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0122_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0141_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-0148_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1230339_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1282934_2808';
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1335201_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1497341_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1499751_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1499752_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1499753_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1500866_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-1516244_2808';                                   
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-251451_2808';                                    
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-542528_2808';                                    
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-545096_2808';                                    
ii=ii+1;ru{ii}='I:\Adela_2808_para Leif\HAM-649829_2808';

resultsFolder='H:\RyRWonLeifs\NEW_RESULTS';

tags={'CON','ISO','CGS','H89','KN93','ZN'};
parrafada={};C=0;
C=C+1;parrafada{C}=['folderName,RR3D,stacksUsed,RR2D,std,stacksUsed,'];
for ii=1:length(ru)
    S=dir(ru{ii});
    
    
    for jj=3:length(S)
        if(S(jj).isdir==1)
            for kk=1:length(tags)
                if(~isempty(strfind(S(jj).name,tags{kk})))
                    try 
%                         processVolume2([ru{ii} '\' S(jj).name],resultsFolder,'l','ch00');
%                         processVolume2([ru{ii} '\' S(jj).name],resultsFolder,'l','ch01');
%                         try
%                             validate2Dwith3D([ru{ii} '\' S(jj).name],resultsFolder,.75,'ch00');
%                             validate2Dwith3D([ru{ii} '\' S(jj).name],resultsFolder,.75,'ch01');
%                         end
                    end
                    try
                       % disp([num2str(jj) '..' num2str(kk)])
                        [r3 m]=overlapping3D([ru{ii} '\' S(jj).name],resultsFolder,1.5,1,{'ch00','ch01'});
                        
                        [r2 s2 n]=overlapping2D([ru{ii} '\' S(jj).name],resultsFolder,1.5,1,{'ch00','ch01'});
                        
                        save([resultsFolder '\' S(jj).name '\ratios.mat'],'r3','r2','s2','m','n');
                        
                        C=C+1;parrafada{C}=[S(jj).name ',' num2str(r3) ',' num2str(m) ',' num2str(r2) ',' num2str(s2) ',' num2str(n) ','];
                    end
                    break;
                end
            end
            
        end
    end
    
end

 muntaCSV([resultsFolder '/allRedRatios.csv'],',',parrafada');
end


