function [ hash_val] = hash_func(random_projection,str,col,base)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

          sum = 0;
          for m = 1:col
                w(m) = base^(random_projection(m)-1);
                if str(random_projection(m)) == 'A'  %1
                    sum = sum + (1 *w(m));
               end
               if str(random_projection(m)) == 'C' %2
                    sum = sum + (2 *w(m));
               end
               if str(random_projection(m)) == 'G' %3
                    sum = sum + (3 *w(m));
               end
               if str(random_projection(m)) == 'T' %4
                    sum = sum + (4 *w(m));
               end
          end
             hash_val = mod(sum,5^base);
end

