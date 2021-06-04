%% 生成图4代码：考虑多种因素的贝叶斯模型预测的疫情趋势图 
%% 包括一下几点内容：
% 1 模型原理
% 2 应用：为什么病毒的传播永远不是指数型的？需要考虑人口空间接触结构
%         可以预测群体免疫阈值（我们提出对于新冠肺炎而言，通过自然免疫实现的群体免疫可能阈值很高，远远高于70%）
%         可以预测人口的空间结构预测疫情的走势（举个例子，均匀分布的结构单一城市群结构和两个城市群之间的结构，其中我们可以给出不同时间感染者空间位置的密度分布！！！用classic5实现）
%         可以预测病毒变异对于疫情的影响（一种可能的情况是群体免疫可能永远无法达到，这反应出病毒与宿主相互进化的过程，可能更接近于真实情况，而且通过注射疫苗的手段甚至是自然群体免疫的手段永远无法根除病毒，用classic4实现）
%% new codes for spreading covid-19 考虑病毒突变后的结果
clc
clear

%% generate correct interaction_matrix 
%% 考虑无症状感染人
clc
clear
for ii = 1
ii
%% 空间分布4个不同的城市，其中ff表示城市1，gg表示城市2，ee表示城市3，hh表示城市4，每个城市2500人，随机分布在100*100的空间区域内，总体种群数量10000人；
ff = 100*rand(2500,2);
gg = 100*rand(2500,2)+110;
ee = 100*rand(2500,2)+220;
hh = 100*rand(2500,2)+330;
dd = [ff;gg;ee;hh];

age = 80*rand(10000,1);%% 考虑到了种群的年龄因素，10000个人年龄从0-80随机分布

distance = zeros(10000,10000);
interaction_M = zeros(10000,10000);
for i = 1: 10000
    for j  = 1:10000
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
for i = 1: 10000
    for j  = 1:10000
        if (distance(i,j) == 0) || (distance(i,j) > 25) %% 使用有距离约束的人口接触矩阵，距离在5以外的点我们认为没有接触概率
            interaction_M(i,j) = 0;
        
        else


              interaction_M(i,j) = min(0.8,10/distance(i,j)^2);
        end

    
    end
end


% 四个城市之间的接触关系
interaction_M(2500,2501) = 0.8;%% 由第一个城市的2500号人和第二个城市的2501号人建立了双向接触关系
interaction_M(2501,2500) = 0.8;
interaction_M(5000,5001) = 0.8;%% 由第二个城市的5000号人和第三个城市的5001号人建立了双向接触关系
interaction_M(5001,5000) = 0.8;
interaction_M(7500,7501) = 0.8;%% 由第三个城市的7500号人和第四个城市的7501号人建立了双向接触关系
interaction_M(7501,7500) = 0.8;
interaction_M(10000,1) = 0.8;%% 由第四个城市的10000号人和第一个城市的1号人建立了单向接触关系





new_matrix = zeros(10000,400);%% 5000代表人数，400代表时间周期
new_matrix(1,1) = 1;%% 一号病人




reverse = new_matrix';
temp_matrix = reverse;
for i = 2:400
    %% 突变矩阵的产生
      mutation_m = zeros(i-1,1);
      for zz = 1:i-1
           mutation_m(zz) = (1-0.005*ii)^(i-zz);%% 认为由病毒突变和抗体衰减因素综合起来的矩阵
      end
      %% 考虑感染发生与病毒的初始入侵概率/数量的关系
      for iii = 1:10000
          temp_matrix(i-1,iii) = reverse(i-1,iii)* 0.1^((1-reverse(i-1,iii))/5); %% 并非所有的接触都会引发感染的发生，它们之间存在非线性关系
      end
      for j = 1:10000



       reverse(i,j) = min((1-sum(reverse(1:i-1,j).*mutation_m)), (1-1/(1+exp(age(j)/25)))^2*(1-sum(reverse(1:i-1,j).*mutation_m))*sum(temp_matrix(i-1,:).*interaction_M(:,j)'));%% 使用考虑多种因素的复杂贝叶斯模型预测感染概率
      end

end




xx(ii) = mean(sum(interaction_M));
original(ii) = max(0,(1-1/xx(ii)));

final_matrix = reverse';
% xx = [1:400];
% plot(xx,sum(final_matrix));
yy(ii) = sum(sum(final_matrix))/10000;

end
%%
new_R_matrix = zeros(10000,10000);
for i = 1:10000
    for j = 1:10000
        new_R_matrix(i,j) = interaction_M(i,j)* 0.1^((1-interaction_M(i,j))/5)*(1-1/(1+exp(age(j)/25)))^2;
    end
end
new_R0 = mean(sum(new_R_matrix));

infection_matrix = temp_matrix';

time_p = [1:400];
xxx = sum(final_matrix);
yyy = sum(infection_matrix);
plot(time_p,xxx,'r');%% 表示每个周期新增抗体阳性人数的变化；
hold on
plot(time_p,yyy,'b');%% 表示每个周期新增核酸检测阳性人数的变化；

%% 从图像中可以看出210世代的时候达到群体免疫水平，下面的代码是计算群体免疫阈值，210是根据上面图来判断的，由于10000个人是随机分布的，每次模拟结果可能会有不同，根据模拟结果选择这个数值。
% 实际产生抗体的比例人数：
antibody_percent = sum(xxx(1:210))/10000;
% 具有传染性的人数，核酸检测阳性比例：
nucl_positive_percent = sum(yyy(1:210))/10000;

% 在这个过程中重复感染的人数以抗体感染判断：
reinfection_rate_calculation


%% 计算平均年龄的变化：
average_age_change_th_time



                