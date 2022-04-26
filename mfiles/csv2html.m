function csv2html(pathtofile)
%pathtofile='F:\Detailed ANALYSIS\Aug27-13 Sparks Linescan-good\inprocess\AllExpsData.csv';


pumpum=fopen(pathtofile,'rt');
A=fread(pumpum,inf,'uchar');
fclose(pumpum);

enters=find(A==10);
if(length(A)~=enters(end))
   enters=[enters; length(A)]; 
end
enters=[1; enters];

sep={',',';'};
AA=char(A)';
[aux,ind]=max([numel(find(AA==sep{1})),numel(find(AA==sep{2}))]);
s=sep{ind};

contingut=cell(1);
cc=0;
cc=cc+1;contingut{cc}='<!DOCTYPE html>';
cc=cc+1;contingut{cc}='<html>';
  cc=cc+1;contingut{cc}='<HEAD>';
  cc=cc+1;contingut{cc}='<style>';
  cc=cc+1;contingut{cc}='table,th,td';
  cc=cc+1;contingut{cc}='{';
  cc=cc+1;contingut{cc}='border:1px solid black;';
  cc=cc+1;contingut{cc}='border-collapse:collapse;';
  cc=cc+1;contingut{cc}='}';
   cc=cc+1;contingut{cc}='</style>';
  cc=cc+1;contingut{cc}='</HEAD>';
  cc=cc+1;contingut{cc}='<BODY>'; 
  cc=cc+1;contingut{cc}='<table>';
  
for ii=2:length(enters)
    
    fra=char(A(enters(ii-1):enters(ii)-1))';
    sec=find(fra==s);
    if(~isempty(sec))
    if(sec(end)~=length(fra)),sec=[sec length(fra)];end
    sec=[1 sec];
     
    frase='<tr>';
    for jj=2:length(sec)
    
        cont=fra(sec(jj-1)+1:sec(jj)-1);
       frase=[frase '<td>' cont '</td>'];
       dad{ii-1,jj-1}=cont;
    end
    frase=[frase '</tr>'];
    cc=cc+1;contingut{cc}=frase;  
    
    end
end

 cc=cc+1;contingut{cc}='</table>';
cc=cc+1;contingut{cc}='</BODY>';
cc=cc+1;contingut{cc}='</html>';


pumpum=fopen([pathtofile(1:end-4) '.html'],'wt');
for ii=1:length(contingut)
   
    fprintf(pumpum,'%s\n',contingut{ii}); 
    
end
fclose(pumpum);


end
