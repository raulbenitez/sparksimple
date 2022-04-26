function [] = comparaExps(DAD,ruta,tags)
FS=18;
N=length(DAD);
%list='k5cFA.geyorz';
 list='rAoTygc5eF.nkz';
%inds=randperm(length(list));
% cols=superjet(N,list(inds(1:N)));
if(N>length(list)),superjet(N,'lines');else cols=superjet(N,list(1:N));end
% properties in the order they come
props={'d2m(um)','d2s(um)','AMP(dF/Fo)','amp(dF/Fo)','Bl(dF/Fo)','RoR(dF/Fo/ms)','mass(ms*dF/Fo)','t2p(ms)','FDHM(ms)','FWHM(um)','tau(ms)'};
startingindex=6;% previous indexes are irrelevant
header='';
for ii=1:length(props)
header=[header ',' props{ii} ',' ];
end
header=[header ','];
%%% 
compDad=zeros(N,2*(1+length(props)));
for ii=1:N
    if(isempty(DAD{ii})),error('Stopped comparison due to empty data.');end
    for jj=1:length(props)
    eval(['v' num2str(jj) '{ii}=DAD{ii}(:,' num2str(jj+startingindex-1) ');']); 
    eval(['compDad(ii,2*jj:2*jj+1)=[mean(v' num2str(jj) '{ii}),std(v' num2str(jj) '{ii})];']);
    end
end
gg=num2cell(compDad);
for ii=1:N
    gg{ii,1}= tags{ii};
end
muntaCSV([ruta '/comparison.csv'],',',header,gg);
 %%%
 scr=get(0,'screensize');
[fil,col]=squareDistrib2(length(props),scr(4)/scr(3));
figure(66);clf;set(gcf,'Menu','none','position',scr,'color','w');
lesX=[];
for ii=1:N
    tags{ii}=[tags{ii} '{\fontsize{' num2str(FS-6) '}{}(' num2str(numel(v1{ii})) ')}'];
lesX=[lesX,ii*ones(1,numel(v1{ii}))];
end
for ii=1:length(props)
    D=[];
    for jj=1:N
eval(['D=[D v' num2str(ii) '{' num2str(jj) '}''];']);
eval(['dists=v' num2str(ii) ';']);
    end
    subNM(fil,col,ii,[.06 .04 .01 .02]);
    fancyBoxplot(D,lesX,cols);
    %gaussPlot(D,lesX,cols);
    statLines(dists);
    grid on;
    ax=gca;
    set(ax,'xtick',[1:N],'xticklabel',tags,'fontname','verdana','fontsize',FS,'box','on','XMinorGrid','on','YMinorGrid','on');ylabel(props{ii});
    set(ax.YLabel,'FontWeight','bold');
end

%saveWysiwyg(66,[ruta '/comparison2.png']);
saveas(66,[ruta '/comparison.png']);
close(66);


end



function [] = statLines(var)
FS=22;
yl=ylim();
N=length(var);
oc0=yl(2);
imp=1;
if(imp==1),pd=-.03;fd=.04;ocF=1.5;ma=.2;else pd=-.06;fd=.05;ocF=2;ma=.5;end
d2=oc0*fd;d3=.05;
oc0=oc0+d2;or=oc0;
ocu=zeros(N,max(0,nchoosek(N,2)+1));
oo=0;
ts=1;
lw=2;
pv='*';
for ii=1:N-1
    %FA=labs{ii};
    for jj=ii+1:N
% FAA=labs{jj};
 if(ts==1),[H,p]=ttest2(var{ii},var{jj});end
% if(ts==2),[H,p]=kstest2(var{ii},var{jj});end
% if(ts==3),p=ranksum(var{ii},var{jj});end
%if(ts==4),p=kruskalwallis(var{ii},var{jj});end
%         disp([FA ' vs ' FAA]);
%         disp(['p-value: ' num2str(p)]);
        if(p<.05)
            if(pv=='*')
           sen=' *';
           if(p<.001),sen=[sen '*'];end
           %if(p<.0001),sen=[sen '*'];end
           %if(p<.00001),sen=[sen '**'];end
            else
               sen=[' p=' niceNumbers(p)]; 
               szc=-5;
            end
           espai=0;or=oc0;
           while(espai==0)
               if(sum(ocu(ii:jj,1+round((or-oc0)/(ocF*d2))))<=1)
                   espai=1;
               else
                   or=or+ocF*d2;
               end
           end
           if(oo==1)
           h=line([or or+d2],[ii+d3 ii+d3]);set(h,'color','k','linewidth',lw);
           h=line([or or+d2],[jj-d3 jj-d3]);set(h,'color','k','linewidth',lw);
           h=line([or+d2 or+d2],[ii+d3 jj-d3]);set(h,'color','k','linewidth',lw);
           h=text(or+d2,pd+(ii+jj)*.5,sen);set(h,'fontsize',FS,'verticalalignment','middle','horizontalalignment','left');
           else
           h=line([ii+d3 ii+d3],[or or+d2]);set(h,'color','k','linewidth',lw);
           h=line([jj-d3 jj-d3],[or or+d2]);set(h,'color','k','linewidth',lw);
           h=line([ii+d3 jj-d3],[or+d2 or+d2]);set(h,'color','k','linewidth',lw);
           h=text((ii+jj)*.5,or+d2*1.5,sen);set(h,'fontsize',FS,'horizontalalignment','center');
           end
           ocu(ii:jj,1+round((or-oc0)/(ocF*d2)))=1;
           
        end
    end
end
end
