%% 用来生成不同时间疫情的空间分布概率图 图4B，可以生成 连续的时间变化



%% 用灰度模型产生，颜色越黑表示感染的概率越大，表示疫情越严重
for i = 1:200
    i
   s = scatter(dd(:,1),dd(:,2),[],1-final_matrix(:,i),'*');
   colormap gray
   caxis([0,1]);
   pause(0.1);
end
% 
% %% 用热度图表示疫情的严重程度 
% 
% for i = 1:200
%     i
%    s = scatter(dd(:,1),dd(:,2),[],1-final_matrix_new(:,i),'*');
%    colormap hsv
%    
%    pause(0.1);
% end
