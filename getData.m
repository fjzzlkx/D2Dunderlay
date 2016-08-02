%直接无脑方式将data里面的cell转换成数组
%资源池数_cellularlink数_D2Dlink数量
a_50_10_100 = cell2mat(data(1,1))
a_50_20_100 = cell2mat(data(1,2))
a_50_30_100 = cell2mat(data(1,3))
a_50_40_100 = cell2mat(data(1,4))

a_50_10_300 = cell2mat(data(1,5))
a_50_20_300 = cell2mat(data(1,6))
a_50_30_300 = cell2mat(data(1,7))
a_50_40_300 = cell2mat(data(1,8))

a_50_10_500 = cell2mat(data(1,9))
a_50_20_500 = cell2mat(data(1,10))
a_50_30_500 = cell2mat(data(1,11))
a_50_40_500 = cell2mat(data(1,12))

%%%%-------------------------%%%%%
a_100_20_100 = cell2mat(data(2,1))
a_100_40_100 = cell2mat(data(2,2))
a_100_60_100 = cell2mat(data(2,3))
a_100_80_100 = cell2mat(data(2,4))

a_100_20_300 = cell2mat(data(2,5))
a_100_40_300 = cell2mat(data(2,6))
a_100_60_300 = cell2mat(data(2,7))
a_100_80_300 = cell2mat(data(2,8))

a_100_20_500 = cell2mat(data(2,9))
a_100_40_500 = cell2mat(data(2,10))
a_100_60_500 = cell2mat(data(2,11))
a_100_80_500 = cell2mat(data(2,12))



%%%%-------------------------%%%%%
a_150_30_100 = cell2mat(data(3,1))
a_150_60_100 = cell2mat(data(3,2))
a_150_90_100 = cell2mat(data(3,3))
a_150_120_100 = cell2mat(data(3,4))
a_150_30_300 = cell2mat(data(3,5))
a_150_60_300 = cell2mat(data(3,6))
a_150_90_300 = cell2mat(data(3,7))
a_150_120_300 = cell2mat(data(3,8))
a_150_30_500 = cell2mat(data(3,9))
a_150_60_500 = cell2mat(data(3,10))
a_150_90_500 = cell2mat(data(3,11))
a_150_120_500 = cell2mat(data(3,12))

a_50_10_100_fix = [a_50_10_100 zeros(MaxTimes1-length(a_50_10_100),1)]; 
a_50_20_100_fix = [a_50_20_100 zeros(MaxTimes1-length(a_50_20_100),1)]; 
a_50_30_100_fix = [a_50_30_100 zeros(MaxTimes1-length(a_50_30_100),1)]; 
a_50_40_100_fix = [a_50_40_100 zeros(MaxTimes1-length(a_50_40_100),1)]; 

a_50_10_300_fix = [a_50_10_300 zeros(MaxTimes1-length(a_50_10_300),1)]; 
a_50_20_300_fix = [a_50_20_300 zeros(MaxTimes1-length(a_50_20_300),1)]; 
a_50_30_300_fix = [a_50_30_300 zeros(MaxTimes1-length(a_50_30_300),1)];  
a_50_40_300_fix = [a_50_40_300 zeros(MaxTimes1-length(a_50_40_300),1)]; 

a_50_10_500_fix = [a_50_10_500 zeros(MaxTimes1-length(a_50_10_500),1)]; 
a_50_20_500_fix = [a_50_20_500 zeros(MaxTimes1-length(a_50_20_500),1)];
a_50_30_500_fix = [a_50_30_500 zeros(MaxTimes1-length(a_50_30_500),1)];
a_50_40_500_fix = [a_50_40_500 zeros(MaxTimes1-length(a_50_40_500),1)];
figure(1)
plot(interationtimes,a_50_10_100_fix,'-o');  hold on;  xlabel('迭代的次数')，ylabel('需要重播的比例');
plot(interationtimes,a_50_20_100_fix,':.o');
plot(interationtimes,a_50_30_100_fix,':o');
plot(interationtimes,a_50_10_300_fix,'-*');
plot(interationtimes,a_50_20_300_fix,':.*');
plot(interationtimes,a_50_30_300_fix,':*');
plot(interationtimes,a_50_10_500_fix,'-+');  
plot(interationtimes,a_50_20_500_fix,':.+');
plot(interationtimes,a_50_30_500_fix,':+');
