function [encode_motif] = encode_motif(str,l_mer)  

s = '';
    for i = 1 : l_mer 
        if str(i) == 'A'  
            s = strcat(s,'00');
        end
        if str(i) == 'C'
            s = strcat(s,'01');
        end
        if str(i) == 'G'
            s = strcat(s,'10');
        end
        if str(i) == 'T'
            s = strcat(s,'11');
        end
    end
encode_motif = s;
end

