function [mutate_offspring] = get_mutation(intermediate_population,mutation_rate,pop_size,l_mer)
%A = 00
%C = 01
%G = 10
%T = 11
%%%%%%%%%%%%%  0.001*100 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
loop_end = (mutation_rate * pop_size);
counter = 1;
mutate_offspring = '';
for j = 1:loop_end
    individual_index   =   randi([1 pop_size],1,1);
    individual_bit_loc =   randi([1 l_mer*2],1,1);
    encod_parent1      =  encode_motif(intermediate_population(individual_index,:),l_mer);
    if encod_parent1(individual_bit_loc) == '1'
        encod_parent1(individual_bit_loc) = '0';
    else
        encod_parent1(individual_bit_loc) = '1';
    end
    mutate_offspring(counter,:) = decode_motif(encod_parent1,l_mer);
end
end

