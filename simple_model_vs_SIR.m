%% 图一生成代码：使用三种不同模型预测的种群感染趋势
%% 包括一下几点内容：
% 1 和SIR模型的比对
% 2 简单模型，不包括复杂地理因素和突变因素

clc
clear



clc
clear
for ii = 1
ii

% generate correct interaction_matrix 
dd = 250*rand(10000,2);



% dd = 250*rand(10000,2);
distance = zeros(10000,10000);
interaction_M = zeros(10000,10000);
for i = 1: 10000
    for j  = 1:10000
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
%% 使用一种无距离约束条件的接触矩阵
for i = 1: 10000
    for j  = 1:10000
        if distance(i,j) == 0  
            interaction_M(i,j) = 0;
        
        else

            interaction_M(i,j) = min(0.8,5/distance(i,j)^2); %%使用一种无距离约束条件的接触矩阵
        end

    
    end
end
%% 使用一种有距离约束条件的接触矩阵
interaction_M_new = zeros(10000,10000);
for i = 1: 10000
    for j  = 1:10000
        if distance(i,j) == 0  || distance(i,j) >=25   %%使用一种有距离约束条件的接触矩阵
            interaction_M_new(i,j) = 0;
        
        else

            interaction_M_new(i,j) = min(0.8,5/distance(i,j)^2);
        end

    
    end
end


new_matrix = zeros(10000,200);%% 10000代表人数，200代表时间周期
new_matrix(1,1) = 1;% 表示疫情开始阶段1号病人出现



reverse = new_matrix';
reverse_new = reverse;

for i = 2:200
    i

      for j = 1:10000




        reverse(i,j) = min(1-sum(reverse(1:i-1,j)),(1-sum(reverse(1:i-1,j)))*sum(reverse(i-1,:).*interaction_M(:,j)'));
        reverse_new(i,j) = min(1-sum(reverse_new(1:i-1,j)),(1-sum(reverse_new(1:i-1,j)))*sum(reverse_new(i-1,:).*interaction_M_new(:,j)'));
      end

end

final_matrix = reverse';
final_matrix_new = reverse_new';

end
%%  使用SIR模型预测的结果代码
%使用的最简单的SIR模型，其中不考虑人口的死亡，认为感染者都可以在一个感染周期后康复，且对病毒完全免疫
yyy = zeros(1,201);
yyy(1) = 1;
sum_yy = 0;
for i = 1:200
    
    R0 = mean(sum(interaction_M));
    sum_yy = sum_yy + yyy(i);
    r = R0*(10000-sum_yy)/10000;
    yyy(i+1) = yyy(i)*r;
 
    
end

time_p = [1:200];
xxx = sum(final_matrix);
xxx_new = sum(final_matrix_new);

plot(time_p,xxx(1:200),'r');% 表示使用无距离约束条件的贝叶斯模型预测的感染曲线
hold on
plot(time_p,xxx_new(1:200),'b');% 表示使用有距离约束条件的贝叶斯模型预测的感染曲线
hold on
plot(time_p,yyy(1:200),'y');% 表示SIR模型预测的感染曲线




                