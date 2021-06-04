overall_anti= 0;

% first 50 generation
for i = 1:10000
    
        if sum(final_matrix(i,1:210)) > 1
            temp_d = 1;
        else
            temp_d = sum(final_matrix(i,1:210));
        end
        overall_anti = overall_anti + temp_d;
end
reinf_rate = 1 - overall_anti/sum(xxx(1:210));
        




