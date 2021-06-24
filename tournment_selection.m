function [ interm_pop ] = tournment_selection(consenses_motif,eval_motif,pop_limit)
% function use for generating intermediate population
for i=1:pop_limit
    r = randi([1 pop_limit],1,2);
    if eval_motif(r(1))>eval_motif(r(2))
        interm_pop(i,:) = consenses_motif(r(1),:);
    else
        interm_pop(i,:) = consenses_motif(r(2),:);
    end
     
end
end

