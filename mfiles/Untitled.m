








 G=fopen('SparkSimple2.m','r');
%cont=textscan(G,'%s');% actual file contents
cont=fread(G,inf,'uint8');% actual file contents
fclose(G);

fins=sort([find(cont==10);find(cont==13)]);
if(fins(end)~=length(cont)),fins=[fins;length(cont)+1];end
prev=1;
parrafada=cell(0);cc=0;c2=0;relevant=[];
for ii=1:length(fins)

    tram=char(cont(prev:fins(ii)-1)');
    try
    if(strcmp(tram(1:12),'%%%%%%%%%%%%'))
        c2=c2+1;
    end
    if(c2==2)
        try
            if(~strcmp(tram(1),'%'))
                relevant=[relevant,cc+1];
            end
        end
    end
    end
    if(~isempty(tram))
    cc=cc+1;
    parrafada{cc}=tram;
    end
    
    prev=fins(ii)+1;
end

for ii=1:length(relevant)
   disp(parrafada{relevant(ii)}); 
end

pumpum=fopen('SparkSimple22.m','wt');
for ii=1:length(parrafada)
   
    fprintf(pumpum,'%s\n',parrafada{ii}); 
    
end
fclose(pumpum);
