function [ evalution_res ] = evaluation_func( c_dna,str,row,n,str_mer,l_mer,increment_str,string_len,pop_limit)
evalution_res = zeros(pop_limit,1);
counter = 1;
h = 0;
min_score = 0;
min = 9999;
for pop_size = 1 : pop_limit
    
    for i=1:row
        increment_str = l_mer;
        for j=1 : ((n - l_mer) + 1)
            counter = 1;
            for l = j:increment_str
                if l <= string_len
                    str(counter) = c_dna(i,l);
                    counter = counter + 1;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for g=1:l_mer
                if str(g) ~= str_mer(pop_size,g)
                    h = h + 1;
                end
            end
            if min > h
                min = h;
            end
            h = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            increment_str = increment_str + 1;
        end
        min_score = min_score + min;
        min = 9999;
    end
    fv = (row * string_len) - min_score;
    evalution_res(pop_size,1) = fv;
    min_score = 0;
end
avg = sum(evalution_res)/pop_limit;
evalution_res(:,1) = evalution_res(:,1)/avg;
[rows cols]   = size(evalution_res);
end


