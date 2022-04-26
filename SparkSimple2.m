function [] = SparkSimple2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SparkSimple2    LASSIElab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%       AlexVL    2013-2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=1;            % Number of experiments to process (set N=0 to analyze all sub-folders in folder)
tagSpk='ch01';  % Relevant channel identifier
RF='ResultsNewest';% Name for results folder

% Detection settings
roiR=1.5;         % ROI radius around spark candidate center (in um)
fs=0;           % Gaussian filter size (in um) for spatial smoothing in images (photon counting issues)
ft=1;           % Temporal filtering (in frames)
dm=2;           % Distance in microns for merging
tm=30;          % Time in miliseconds for merging
nc=200;          % Number of candidates

% Filtering parameters
iT=0.2;         % Spark relative intensity threshold (maximum value)
Bl=1;            % Spark previous baseline (maximum value)
Fd=[10 200];      % Full duration at half maximum in ms (valid interval)
ta=[10 200];     % Spark decay constant tau in ms (valid interval)
t2=[0 200];      % Time to peak for spark uprise in ms (valid interval)
r2=0.4;           % Spark R² on decay fit (minimum value)

% Aditional settings
ss=1;            % Signal sign (1 for positive, -1 for negative)
fb=0;            % Force baseline (if N>1 first folder will force baseline value for the rest)
ys=4.5;          % Y-scale max val for outliers
bf=1;            % Brightness factor for output frames
sr=1;            % Show rejected sparks in film
sp=1;            % Starting pont (1-detection,2-features,3-filter,4-rois,5-export)
expNam={'CON','DRUG1','DRUG2','DRUG3','DRUG4','DRUG5','DRUG6'};% exp names for plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off','all')
%load('ultimpaquet.mat');
addpath('./mfiles');
%if(fb==0),fb=-1;end% per poder reciclar variables
save('./mfiles/RF.mat','RF');gg=ys;
save('./mfiles/ys.mat','ys');
prefix='yData';save('./mfiles/prefix.mat','prefix');
if(N>0)
    folder=' ';try load('ultimpaquet.mat');folder=rutes{1};end
    for ii=1:N
        try folder=uigetdir(folder);catch folder=uigetdir(' ');end
        rutes{ii}=folder;folder=folder(1:find(folder=='\',1,'last')-1);
        if(N>1)
            if(ii==1),ensenya('Adding to processing pipeline:');end
        disp(['        ' num2str(ii) '-. ' rutes{ii}]);
        end
    end
else
    folder=' ';try load('ultimpaquet.mat');folder=rutes{1};end
    try folder=uigetdir(folder);catch folder=uigetdir(' ');end
    s=dir(folder);
    cc=0;
    for ii=3:length(s)
        if(s(ii).isdir)
            cc=cc+1;
            rutes{cc}=[folder '/' s(ii).name];
        end
    end
    N=length(rutes);
end
save('./mfiles/ultimPaquet','rutes');

try gg=urlread('http://leica.upc.es/counters/ss.php'); end
firstfolder=rutes{1};
for ii=1:N
    ensenya(['Starting with folder ' rutes{ii}],superjet(1,'e'));
    try
        folder=rutes{ii};
        S=[dir([folder filesep '*' tagSpk '*.jpg']) dir([folder filesep '*' tagSpk '*.tif'])];
        if isempty(S)
            S=dir([folder filesep '*' tagSpk '*.data']);
        end
        DAT=[];
        if(isempty(S)),error('No images found!');end
        try
        aux=imread([folder filesep S(1).name]);
        disp(['         Physical size: ' niceNums(size(aux,1),0) ' x ' niceNums(size(aux,2),0) ' x ' niceNums(length(S),0) '.']);        
        catch
        end
        if(sp==1)
             try rmdir([folder filesep RF],'s'); end
        end
             Rfol=[folder filesep RF];
            if(exist(Rfol,'dir')~=7),mkdir(Rfol);end
            if(exist([Rfol '/detFilm'],'dir')~=7),mkdir([Rfol '/detFilm']);end
            if(exist([Rfol '/extraFigs'],'dir')~=7),mkdir([Rfol '/extraFigs']);end
            if(sp<2),delete([Rfol '/extraFigs/*.*']);end
            %------------% 1 detection
            if(sp==1)   
               %if((ii==1)),fba=-1;else fba=fb;end% 
               if(ii==1),fo=fb;else,if(fb==0),fo=0;end;end
               fo = sparkDetect(folder,Rfol,S,dm,tm,roiR,fs,nc,ft,ss,fo);
               %if((ii==1)&&(fb~=-1)),fb=fo;end% reciclem el flag de force baseline i assignem el valor a arrossegar
            else
                load([Rfol '/zMetaData.mat']);
            end
            %------------% 2 features
            if(sp<=2)
            sparkFeatures(folder,Rfol,roiR);
            end
            %------------% 3 filter
            if(sp<=3)
            sparkFiltSimple(folder,Rfol,iT,ta,t2,r2,Bl,Fd);            
            end
            %------------% 4 rois
            if(sp<=4)
            clusterROIs(folder,Rfol,S,sr,bf);    
            end
            load([Rfol '/zData2.mat']);if(exist('spkF','var')==0),spkF=SPARKS;end            
            %------------% 5 export
            if(sp<=5)&&(isempty(spkF)==0)
                [dades] = guardaDades(spkF,ROIs,folder,Rfol);
            end
            DAT=[DAT;dades];
            ensenya('Folder finished.');
            
%         end
        DAD{ii}=DAT;
    catch er;ensenya(['Skipping due to error: ' er.message ' In ' er.stack(1).name '.m (line ' num2str(er.stack(1).line) ')'],superjet(1,'r'));end
    disp(' ');
end
try
comparaExps(DAD,[firstfolder,'\',RF],expNam);
end
close all;
end