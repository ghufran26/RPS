 t =10;
 k = 7;
 theta = 4;
 k_pos = [1 2 3 4 5 6 7];
 n=length(k_pos);
 nr=round(n/1.5);
 string_len = 100;
 decode_kmers = [1,2,3,4];
 random_projection = k_pos(randperm(n,nr));
 [row1 col1] = size(random_projection);
 disp(random_projection);
 base = k;
 hashfuc = 0;
 sum = 0;
 w = [];
%AGA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_dna = [
'ACAAAACCCATCGTAGTCCTTAGACTTGGGACACTTACACCTGCAGCGCGCGCATGTGGAAATAGAGGCCAAGTTCGATCCGTACTCCGACGTACGATGC';
'AACAGTGTGGATGTGACGAGATTCATTTATACCCTTCGCGCGCCGGACTGGCCTCGGCAAGGCGCGGCGGTGAACAAGCAATTGACAACTAACCACCGTG';
'TATTCGTTATGGCATAAGGCAGTTTAAGTCGAGACAATAGGGCTCGCAATACACAGTTTACCGCATATTGCCCTAACTGACAAACTGTGATCGACCACTA';
'GCCATGCCATTGCCTCTTAGATACCGCGATACAGTGATTATGAAAGGTTTGTGGGGCATGGCTACGACTTGTTCAGCTACGTCCGAGGGCAGAAACTTAT';
'CCCCATTTGTATGTTGACCTATCTACTACCGATCCCCGGAGGTTAAGTAGGTTGTGAGATGCGGGAGAGGTTCTCGATCTTCCCGTGGGACGTCAACCTT';
'TCCCTTGATAAAGCATCCCGCTCGGGTATGGCAGTGAGTACGCCTTCTGAATTGTGCTATCCTTCGTCCTTATCAAAGCTTGCTACCAATAATTAGGATT';
'ATTGCCTTGCGACAGACTTCCTACTCACACTCCCTCACATTGAGCTACTCGATGGGCGATTAGCTTGACCCGCTCTGTAGGGTCGCGACTACGTGAGCTA';
'GGGCTCCGGACTGGGCTGTATAGTCGAGTCTGATCTCGCCCCGACAACTGCAAACCCCAACTTATTTAGATAACATGGTTAGCCGAAGTTGCACGGGGTG';
'CCGACCGTGGACTCCTCCCCGGGTGTGGCTCGTTCATCTGACAACATGCAAGCGCTACCACCATCGATTGATTCAGCGGACGGTGTTGTTGTCATAGATT';
'CGGCACATTTCTCTTGTAGGTGTGAAATCACTTAGGTTCGCGCCGTAGTCTTATGGCAAAACCGATGGACTATGTTTCGGGTAGCACCAGGAGTCTGTAG';
];

%fileread('E:\MSCS-1\Advance Algo\project\testexam.txt');
n = 100;
l_mer = k;
num_B = zeros(1,5^base);
[row,col] = size(c_dna);

increment_str = k;
%%%%%%%%%%%Bucket array %%%%%%%

   Bucket_array =cell(5^base,n); 
   bucket_index = [];
   buk_counter = 1;
   Main_consenses_String = 'ggg';
   consenses_counter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['a','b','c','d','e'];
counter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start RPS Algorithm %%%%%%%%%%%%%%%%%%%
for i=1:row
    for j=1 : ((n - l_mer) + 1) 
         counter = 1;
          for l = j:increment_str
           if l <= string_len
            str(counter) = c_dna(i,l);
            counter = counter + 1;
           end
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
             hashfuc = hash_func(random_projection,str,col1,5);
             num_B(hashfuc) = num_B(hashfuc) + 1;
             sum = 0;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             increment_str = increment_str + 1;
    end
end

%disp(num_B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

increment_str = k;
for i=1:row
    for j=1 : ((n - l_mer) + 1) 
         counter = 1;
          for l = j:increment_str
              if l <= string_len
                str(counter) = c_dna(i,l);
                counter = counter + 1;
              end
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             hashfuc = hash_func(random_projection,str,col1,5);
             if num_B(hashfuc) >= theta
                  buk = getIndex(Bucket_array,hashfuc);                  
                  Bucket_array{hashfuc}{buk} = str;
                  st = Bucket_array{hashfuc}{buk};
                  buk_counter = buk_counter + 1;
                  bucket_index(buk_counter)=hashfuc;
             end          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          increment_str = increment_str + 1;
    end
end

%%%%%%%% string size is 5 now %%%

unique_bucket = unique(bucket_index);
[row_c col_c] = size(unique_bucket);
 for i = 2:col_c
     
     buk_index = getIndex(Bucket_array,unique_bucket(i));
     buk_index = buk_index - 1;
     count_profile = zeros(4,k);

      for j = 1:buk_index
          str1 = Bucket_array{unique_bucket(i)}{j};  
          [row3 col] = size(str1);
          for y = 1:k
               disp(str1(y));
              if str1(y) == 'A'
                 count_profile(1,y) = count_profile(1,y) + 1; 
               end
               if str1(y) == 'C'
                 count_profile(2,y) = count_profile(2,y) + 1; 

               end
               if str1(y) == 'G'
                   count_profile(3,y) =count_profile(3,y) + 1;
               end
               if str1(y) == 'T'
                   count_profile(4,y) =count_profile(4,y) + 1;
               end
          end
      end
      count_profile = count_profile/buk_index;
      disp('aaaa');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      increment_str = k;
      start_index = 1;
      end_index   = k;
      consenses_string='akak';%(5^base,col);

      for ii=1:row
          start_index = 1;
          end_index   = k;
          increment_str = k;
           max = 0;
           
         for jj=1 : ((n - l_mer) + 1) 
           counter = 1;
          for l = jj:increment_str
              if l <= string_len
                    str2(counter) = c_dna(ii,l);
                    counter = counter + 1;
              end
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          [str_row,str_col] = size(str2);
          prob = 1;
            
          for str_i = 1:k
              if str2(str_i) == 'A'
                  prob = prob * count_profile(1,str_i);   
              end
              if str2(str_i) == 'C'
                 prob = prob * count_profile(2,str_i);
              end
              if str2(str_i) == 'G'
                 
                  prob = prob * count_profile(3,str_i);
              end
              if str2(str_i) == 'T'
                  prob = prob * count_profile(4,str_i);
              end
          end
          if max < prob
             max = prob;
             start_index = jj;
             end_index   = increment_str;
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          increment_str = increment_str + 1;
         end
         str4='';
         start_counter = 1;
         for jjj=start_index:end_index
            if jjj ~= (n + 1)
                consenses_string(ii,start_counter) = c_dna(ii,jjj);
                start_counter = start_counter + 1;
            end
           % disp(c_dna(ii,jjj));
         end
        
      end

      [row_s,col_s] = size(consenses_string);
      final_count_profile = zeros(4,k);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i_row=1:row_s
            for str_j = 1:k
                  if consenses_string(i_row,str_j) == 'A'
                      final_count_profile(1,str_j)  = final_count_profile(1,str_j) + 1;  
                  end
                  if consenses_string(i_row,str_j) == 'C'
                      final_count_profile(2,str_j)  = final_count_profile(2,str_j) + 1;  
                  end
                  if consenses_string(i_row,str_j) == 'G'
                      final_count_profile(3,str_j)  = final_count_profile(3,str_j) + 1; 
                  end
                  if consenses_string(i_row,str_j) == 'T'
                      final_count_profile(4,str_j)  = final_count_profile(4,str_j) + 1;
                  end
            end
        end
       max = 0;
       index_gg = [];
       for hh = 1:k
           index_gg = 1;
           max = 0;
           for gg = 1:4
               if max < final_count_profile(gg,hh)
                  max = final_count_profile(gg,hh);
                  index_gg = gg;
               end
               
           end
           if index_gg == 1
             Main_consenses_String(consenses_counter,hh) = 'A';  
           end
           if index_gg == 2
             Main_consenses_String(consenses_counter,hh) = 'C';  
           end
           if index_gg == 3
             Main_consenses_String(consenses_counter,hh) = 'G';  
           end
           if index_gg == 4
             Main_consenses_String(consenses_counter,hh) = 'T';  
           end
       end
       consenses_counter = consenses_counter + 1; 
     % disp(m);
     % disp(l);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end
      disp(consenses_counter);
 %%%%%%%%%%%%%%%%%%%%%%%%%%consenses string %%%%%%%%%%%%%%%
 disp('Consenses String generated from RPS');
 disp(Main_consenses_String);
 %{
 %%%%%%%%%%%%%%%%%%%%%%%%  End Algorithm RPS %%%%%%%%%%%%%%%%%%%%%%%%%
%}
 
 

