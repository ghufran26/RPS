function [crossover_offspring] = get_crossover(intermediate_population,crossover_rate,pop_size,l_mer)
%A = 00
%C = 01
%G = 10
%G = 11
%%%%%%%%%%%%%  0.30*100 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
loop_end = (crossover_rate * pop_size)/2;
counter = 1;
for j = 1:loop_end
    r = randi([1 pop_size],1,2);
    parent1(:) = intermediate_population(r(1),:);
    parent2(:) = intermediate_population(r(2),:);
    encod_parent1 =  encode_motif(parent1,l_mer);
    encod_parent2 =  encode_motif(parent2,l_mer);
    center  =(l_mer*2)/2;
    offspring1 = '';
    offspring2 = '';
    
    for i=1:(l_mer*2)
        if i <= center
            offspring1 = strcat(offspring1,encod_parent1(i));
            offspring2 = strcat(offspring2,encod_parent2(i));
        else
            offspring1 = strcat(offspring1,encod_parent2(i));
            offspring2 = strcat(offspring2,encod_parent1(i));
        end
    end
    decode_offspring1 =  decode_motif(offspring1,l_mer);
    decode_offspring2 =  decode_motif(offspring2,l_mer);
    
    crossover_offspring(counter,:) = decode_offspring1(:);
    counter = counter + 1;
    crossover_offspring(counter,:) = decode_offspring2(:);
    counter = counter + 1;
end
end

