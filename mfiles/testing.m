

gg=mean(abs(diff(volum,1,3)),3);
ed=linspace(0,max(max(gg)),100);
h=histc(reshape(gg,[1,numel(gg)]),ed);
figure,bar(ed,h);hold on;

g=[];b=[];
for ii=1:length(spkF)
    var=max(spkF(ii).timesignal)-min(spkF(ii).timesignal);
    if(spkF(ii).good==1)
       g=[g var];
    else
       b=[b var]; 
    end
    
end
h=histc(g,ed);bar(ed,h);
h=histc(b,ed);bar(ed,h);