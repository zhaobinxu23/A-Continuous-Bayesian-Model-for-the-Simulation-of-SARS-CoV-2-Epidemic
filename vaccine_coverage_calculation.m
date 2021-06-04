%% 图三生成代码：使用简单的贝叶斯模型预测最佳疫苗接种比例

% 使用简单模型计算疫苗接种的最佳覆盖率；不考虑突变效应

clc
clear

%% generate correct interaction_matrix

clc
clear
for ii = 1:10
ii 
%% ii表示疫苗接种的比例，其比例= 1-ii/10；即90%，80%，70%.......0%十个不同的接种梯度，造成在250*250面积下实际非接种人口的减少，通过这种减少会影响R0
dd = 250*rand(1000*ii,2);



% dd = 250*rand(1000*ii,2);
distance = zeros(1000*ii,1000*ii);
interaction_M = zeros(1000*ii,1000*ii);
for i = 1: 1000*ii
    for j  = 1:1000*ii
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
for i = 1: 1000*ii
    for j  = 1:1000*ii
        if distance(i,j) == 0  
            interaction_M(i,j) = 0;
        
        else

            interaction_M(i,j) = min(0.6,15/distance(i,j)^2);% 使用无距离限制的种群接触密度矩阵
        end

    
    end
end


new_matrix = zeros(1000*ii,100);%% 1000*ii代表人数，100代表时间周期
new_matrix(1,1) = 1;% 表示疫情开始阶段1号病人出现



reverse = new_matrix';


for i = 2:100
    i

      for j = 1:1000*ii




        reverse(i,j) = min(1-sum(reverse(1:i-1,j)),(1-sum(reverse(1:i-1,j)))*sum(reverse(i-1,:).*interaction_M(:,j)'));%% 使用贝叶斯方法确定每个个体的具体时间的感染概率
        
      end

end

final_matrix = reverse';
xx(ii) = mean(sum(interaction_M)); %% xx表示不同接种率造成的R0值的变化，理论上是线性相关，因为我们的模型使用了随机人口分布，因此结果没有呈现出完全的线性相关性；
original(ii) = max(0,(1-1/xx(ii)));%% 表示根据R0值推导的群体免疫阈值


yy(ii) = sum(sum(final_matrix))/(1000*ii);%% yy表示未接种疫苗人群感染的概率
zz(ii) = (sum(sum(final_matrix)) + (10000-1000*ii))/10000;%% zz表示要实现群体免疫最终人群的抗体阳性率

end





time = [1:10];
plot(time,xx,'r');
hold on
plot(time,original,'b');
hold on
plot(time,yy,'y');
hold on
plot(time,zz,'c');






                