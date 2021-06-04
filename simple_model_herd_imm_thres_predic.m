%% 图2生产代码：使用不同方法预测出的群体免疫阈值比较



%  使用简单模型计算不同R0值的病毒群体免疫的阈值；不考虑突变效应


clc
clear

%% generate correct interaction_matrix 

clc
clear
for ii = 1:10
    %% ii 表示不同的种群接触频率，对应于不同的R0值，我们使用了250*250范围内随机分布10000个个体
ii
dd = 250*rand(10000,2);



% dd = 250*rand(10000,2);
distance = zeros(10000,10000);
interaction_M = zeros(10000,10000);
for i = 1: 10000
    for j  = 1:10000
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
for i = 1: 10000
    for j  = 1:10000
        if distance(i,j) == 0  
            interaction_M(i,j) = 0;
        
        else

            interaction_M(i,j) = min(0.6,3*ii/distance(i,j)^2); %% 采用了无距离约束条件的贝叶斯模型
        end

    
    end
end


new_matrix = zeros(10000,100);%% 10000代表人数，100代表时间周期
new_matrix(1,1) = 1;% 表示疫情开始阶段1号病人出现





reverse = new_matrix';


for i = 2:100
    i

      for j = 1:10000




        reverse(i,j) = min(1-sum(reverse(1:i-1,j)),(1-sum(reverse(1:i-1,j)))*sum(reverse(i-1,:).*interaction_M(:,j)'));
        
      end

end

final_matrix = reverse';
xx(ii) = mean(sum(interaction_M));%% xx（ii）表示根据人口接触频率计算出的R0值
original(ii) = max(0,(1-1/xx(ii)));%% original(ii)表示根据R0值推导出的群体免疫阈值  关系 ：群体免疫阈值 = 1-1/R0；


yy(ii) = sum(sum(final_matrix))/10000;

end



%%  使用SIR模型预测的群体免疫阈值
%使用的最简单的SIR模型，其中不考虑人口的死亡，认为感染者都可以在一个感染周期后康复，且对病毒完全免疫

for k = 1:10
    yyy = zeros(1,101);
    yyy(1) = 1;
    sum_yy = 0;
    R0 = xx(k);% 对应于不同病毒的R0值
for i = 1:100
    
    
    sum_yy = sum_yy + yyy(i);
    r = R0*(10000-sum_yy)/10000;

    yyy(i+1) = min(yyy(i)*r,(10000-sum_yy));
end
 yy_sir(k) = sum_yy/10000;% yy_sir（k）表示在根据SIR模型计算，最终抗体阳性的比例
    
end


time = [1:10];
plot(time,xx,'r');%% 对应于不同状况下的R0值
hold on
plot(time,original,'b');%% 对应于使用简单公式根据R0推导出的群体免疫阈值
hold on
plot(time,yy,'r');%% 对应于使用贝叶斯模型计算出的群体免疫阈值
hold on
plot(time,yy_sir,'c'); % yy_sir表示在根据SIR模型计算，最终抗体阳性的比例，也就是群体免疫的阈值



                