for i = 1:400
    average_new = 0;
    for j = 1:10000
        average_new = average_new + age(j)*final_matrix(j,i);
    end
    aver_age(i) = average_new/xxx(i);
end
time_l = [1:400];
plot(time_l,aver_age,'b')

        