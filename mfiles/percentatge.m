function [] = percentatge(ara,tots)

    

   numeret=round(100*ara/tots);
   n0=numeret;
   numeret=num2str(numeret);
   if(n0<100),numeret=['0' numeret];end
   if(n0<10),numeret=['0' numeret];end
    fprintf(['\b\b\b\b\b ' numeret '%%']);  
   
   
   if ara==tots
       %fprintf(1,'\b\b\b\b\b\n');
       fprintf(['\b\b\b\b']);fprintf(['\n']);
   end
end

