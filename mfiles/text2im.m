function [frase] = text2im(txt)

load('charset.mat');

frase=[];
for ii=1:length(txt)
    ind=strfind(index,txt(ii));
    if(~isempty(ind))
   frase=[frase chars{ind(1)}]; 
    end
end


end