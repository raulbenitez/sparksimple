function [dades] = csv2cell(ruta)



          G=fopen(ruta,'rt');
          cont=(fread(G,inf,'uchar'))';% l
          fclose(G);
 if(numel(find(cont==44))>numel(find(cont==59))),
     separador=44;
 else
     separador=59;
 end
%separador=44;% comma
%separador=59;% semicolon        
          
          enters=sort([find(cont==10) find(cont==13)]);

          dades=cell(numel(enters),1);
          prev=1;
          for ii=1:numel(enters)
              tram=cont(prev:enters(ii)-1);
              seps=find(tram==separador);
              pr2=1;
              
              for jj=1:numel(seps)
                  
                 dades{ii,jj}=char(tram(pr2:seps(jj)-1));
                 pr2=seps(jj)+1;
              end
              %dades{ii,1}=tram;
              if(pr2<length(tram)),
                  if(numel(seps)==0),jj=0;end
                dades{ii,jj+1}=char(tram(pr2:end)); 
              end
              prev=enters(ii)+1;
          end
          


end