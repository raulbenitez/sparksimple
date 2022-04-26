function [] = pintaDades(DATA)






figure;set(gcf,'position',[534   673   683   649],'color','w');
%%%%  i  X Y t AMP amp bl ror t2p FDHM FWHM tau
noms={'Spark Amplitude','Full duration at half maximum','Full width at half maximum'};
labs={'(F-Fo)/Fo','ms','um'};
DATA(:,11)=DATA(:,11).*.07;
DATA(:,13)=DATA(:,6)./DATA(:,7);
lloc=[13 10 11];
xlims=[0 5;0 500;0 6];
eds=[100,50,100]+1;
for ii=1:length(noms)
   
    
vec=DATA(:,lloc(ii));
ed=linspace(xlims(ii,1),xlims(ii,2),eds(ii));
h=histc(vec,ed);

subplot(1,3,ii);set(gca,'fontsize',14);

bar(ed,h,'k');title(noms{ii});xlim(xlims(ii,:));
       
%     gg=ylim;gc=gg(2)/1000;
%    h=line([ed(2) ed(2)]/10,gg);set(h,'color','k');
%    gg=xlim;
%    h=line(gg,[gc gc]);set(h,'color','k'); 
  xlabel(labs{ii});
end




% 
% 8 features
% 
% figure;set(gcf,'position',[534   673   683   649],'color','w');
% %%%%  i  X Y t AMP amp bl ror t2p FDHM FWHM tau
% noms={'Absolute amplitude','Relative amplitude','Local Baseline','Rate of Rise','time to peak','Full duration at half maximum','Full width at half maximum','Decay constant'};
% lloc=[5 6 7 8 9 10 11 12];
% xlims=[0 1;0 1;0 1;0 .02;0 500;0 500;0 100;0 500];
% eds=[100,100,100,100,100,50,100,100]+1;
% for ii=1:length(noms)
%    
%     
% vec=DATA(:,lloc(ii));
% ed=linspace(xlims(ii,1),xlims(ii,2),eds(ii));
% h=histc(vec,ed);
% 
% subplot(2,4,ii);
% 
% bar(ed,h,'k');title(noms{ii});xlim(xlims(ii,:));
%        
%     gg=ylim;gc=gg(2)/1000;
%    h=line([ed(2) ed(2)]/10,gg);set(h,'color','k');
%    gg=xlim;
%    h=line(gg,[gc gc]);set(h,'color','k'); 
%   
% end
% 



end