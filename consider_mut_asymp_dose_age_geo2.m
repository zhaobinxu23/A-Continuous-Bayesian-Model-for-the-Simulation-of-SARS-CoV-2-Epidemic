%% ����ͼ4���룺���Ƕ������صı�Ҷ˹ģ��Ԥ�����������ͼ 
%% ����һ�¼������ݣ�
% 1 ģ��ԭ��
% 2 Ӧ�ã�Ϊʲô�����Ĵ�����Զ����ָ���͵ģ���Ҫ�����˿ڿռ�Ӵ��ṹ
%         ����Ԥ��Ⱥ��������ֵ��������������¹ڷ��׶��ԣ�ͨ����Ȼ����ʵ�ֵ�Ⱥ�����߿�����ֵ�ܸߣ�ԶԶ����70%��
%         ����Ԥ���˿ڵĿռ�ṹԤ����������ƣ��ٸ����ӣ����ȷֲ��Ľṹ��һ����Ⱥ�ṹ����������Ⱥ֮��Ľṹ���������ǿ��Ը�����ͬʱ���Ⱦ�߿ռ�λ�õ��ܶȷֲ���������classic5ʵ�֣�
%         ����Ԥ�ⲡ��������������Ӱ�죨һ�ֿ��ܵ������Ⱥ�����߿�����Զ�޷��ﵽ���ⷴӦ�������������໥�����Ĺ��̣����ܸ��ӽ�����ʵ���������ͨ��ע��������ֶ���������ȻȺ�����ߵ��ֶ���Զ�޷�������������classic4ʵ�֣�
%% new codes for spreading covid-19 ���ǲ���ͻ���Ľ��
clc
clear

%% generate correct interaction_matrix 
%% ������֢״��Ⱦ��
clc
clear
for ii = 1
ii
%% �ռ�ֲ�4����ͬ�ĳ��У�����ff��ʾ����1��gg��ʾ����2��ee��ʾ����3��hh��ʾ����4��ÿ������2500�ˣ�����ֲ���100*100�Ŀռ������ڣ�������Ⱥ����10000�ˣ�
ff = 100*rand(2500,2);
gg = 100*rand(2500,2)+110;
ee = 100*rand(2500,2)+220;
hh = 100*rand(2500,2)+330;
dd = [ff;gg;ee;hh];

age = 80*rand(10000,1);%% ���ǵ�����Ⱥ���������أ�10000���������0-80����ֲ�

distance = zeros(10000,10000);
interaction_M = zeros(10000,10000);
for i = 1: 10000
    for j  = 1:10000
        distance(i,j) = (dd(i,1) -dd(j,1))^2 + (dd(i,2) -dd(j,2))^2;
    end
end
for i = 1: 10000
    for j  = 1:10000
        if (distance(i,j) == 0) || (distance(i,j) > 25) %% ʹ���о���Լ�����˿ڽӴ����󣬾�����5����ĵ�������Ϊû�нӴ�����
            interaction_M(i,j) = 0;
        
        else


              interaction_M(i,j) = min(0.8,10/distance(i,j)^2);
        end

    
    end
end


% �ĸ�����֮��ĽӴ���ϵ
interaction_M(2500,2501) = 0.8;%% �ɵ�һ�����е�2500���˺͵ڶ������е�2501���˽�����˫��Ӵ���ϵ
interaction_M(2501,2500) = 0.8;
interaction_M(5000,5001) = 0.8;%% �ɵڶ������е�5000���˺͵��������е�5001���˽�����˫��Ӵ���ϵ
interaction_M(5001,5000) = 0.8;
interaction_M(7500,7501) = 0.8;%% �ɵ��������е�7500���˺͵��ĸ����е�7501���˽�����˫��Ӵ���ϵ
interaction_M(7501,7500) = 0.8;
interaction_M(10000,1) = 0.8;%% �ɵ��ĸ����е�10000���˺͵�һ�����е�1���˽����˵���Ӵ���ϵ





new_matrix = zeros(10000,400);%% 5000����������400����ʱ������
new_matrix(1,1) = 1;%% һ�Ų���




reverse = new_matrix';
temp_matrix = reverse;
for i = 2:400
    %% ͻ�����Ĳ���
      mutation_m = zeros(i-1,1);
      for zz = 1:i-1
           mutation_m(zz) = (1-0.005*ii)^(i-zz);%% ��Ϊ�ɲ���ͻ��Ϳ���˥�������ۺ������ľ���
      end
      %% ���Ǹ�Ⱦ�����벡���ĳ�ʼ���ָ���/�����Ĺ�ϵ
      for iii = 1:10000
          temp_matrix(i-1,iii) = reverse(i-1,iii)* 0.1^((1-reverse(i-1,iii))/5); %% �������еĽӴ�����������Ⱦ�ķ���������֮����ڷ����Թ�ϵ
      end
      for j = 1:10000



       reverse(i,j) = min((1-sum(reverse(1:i-1,j).*mutation_m)), (1-1/(1+exp(age(j)/25)))^2*(1-sum(reverse(1:i-1,j).*mutation_m))*sum(temp_matrix(i-1,:).*interaction_M(:,j)'));%% ʹ�ÿ��Ƕ������صĸ��ӱ�Ҷ˹ģ��Ԥ���Ⱦ����
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
plot(time_p,xxx,'r');%% ��ʾÿ�����������������������ı仯��
hold on
plot(time_p,yyy,'b');%% ��ʾÿ�������������������������ı仯��

%% ��ͼ���п��Կ���210������ʱ��ﵽȺ������ˮƽ������Ĵ����Ǽ���Ⱥ��������ֵ��210�Ǹ�������ͼ���жϵģ�����10000����������ֲ��ģ�ÿ��ģ�������ܻ��в�ͬ������ģ����ѡ�������ֵ��
% ʵ�ʲ�������ı���������
antibody_percent = sum(xxx(1:210))/10000;
% ���д�Ⱦ�Ե����������������Ա�����
nucl_positive_percent = sum(yyy(1:210))/10000;

% ������������ظ���Ⱦ�������Կ����Ⱦ�жϣ�
reinfection_rate_calculation


%% ����ƽ������ı仯��
average_age_change_th_time



                