function [decode_motif] = decode_motif(str,l_mer)

s = '';
 
   % 'ATC' 001001
    
    cant_str = '';
    for i = 1 : l_mer * 2
        cant_str = strcat(cant_str,str(i));
        if mod(i,2) == 0
            if cant_str == '00'  
                s = strcat(s,'A');
            end
            if cant_str == '01'
                s = strcat(s,'C');
            end
            if cant_str == '10'
                s = strcat(s,'G');
            end
            if cant_str == '11'
                s = strcat(s,'T');
            end   
            cant_str = '';
        end
    end
    decode_motif = s;
end

