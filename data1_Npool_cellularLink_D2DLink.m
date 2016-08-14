%数据含义 a_资源池数量_cellularlink数量_D2Dlink数量
%迭代次数            1         2         3         4         5         6         7         8         9        10       11        12
a_50_10_100   =   [0.0727    0.0182         0];
a_50_20_100   =   [0.1667    0.0417         0];
a_50_30_100   =   [0.2000    0.0692         0];
a_50_40_100   =   [0.3071    0.1571    0.0500    0.0143         0];
a_50_10_300   =   [0.2531    0.0625    0.0125         0];
a_50_20_300   =   [0.3147    0.1147    0.0618    0.0235         0];
a_50_30_300   =   [0.3694    0.1861    0.1222    0.0639    0.0167    0.0056         0];
a_50_40_300   =   [0.4184    0.1974    0.1579    0.0947    0.0579    0.0184    0.0053         0];
a_50_10_500   =   [0.4075    0.1415    0.0509    0.0094         0];
a_50_20_500   =   [0.5054    0.2232    0.1089    0.0643    0.0375    0.0143         0];
a_50_30_500   =   [0.5407    0.2678    0.1644    0.1339    0.1136    0.0966    0.0695    0.0322    0.0119         0];
a_50_40_500   =   [0.5677    0.3161    0.2032    0.1758    0.1516    0.1242    0.1032    0.0742    0.0581    0.0290    0.0097         0];
a_100_20_100  =   [0.0364         0];
a_100_40_100  =   [0.1417         0];
a_100_60_100  =   [0.2154    0.0231         0];
a_100_80_100  =   [0.2071    0.0857    0.0143         0];
a_100_20_300  =   [0.1500    0.0219         0];
a_100_40_300  =   [0.2029    0.0647    0.0118         0];
a_100_60_300  =   [0.3417    0.1222    0.0583    0.0056         0];
a_100_80_300  =   [0.3605    0.1474    0.0474    0.0105         0];
a_100_20_500  =   [0.2774    0.0585    0.0132         0];
a_100_40_500  =   [0.3929    0.1179    0.0464    0.0107         0];
a_100_60_500  =   [0.4169    0.1678    0.0881    0.0305    0.0068         0];
a_100_80_500  =   [0.4306    0.2129    0.1323    0.0790    0.0306         0];
a_150_30_100  =   [0.0636    0.0182         0];
a_150_60_100  =   [0.1083         0];
a_150_90_100  =   [0.1538    0.0154         0];
a_150_120_100 =   [0.1714    0.0286         0];
a_150_30_300  =   [0.1031    0.0094         0];
a_150_60_300  =   [0.2353    0.0500         0];
a_150_90_300  =   [0.2528    0.0667    0.0167         0];
a_150_120_300 =   [0.3053    0.1053    0.0316    0.0053         0];
a_150_30_500  =   [0.2377    0.0302         0];
a_150_60_500  =   [0.2857    0.0750    0.0125         0];
a_150_90_500  =   [0.3169    0.1220    0.0373    0.0068         0];
a_150_120_500 =   [0.4210    0.1710    0.0919    0.0258         0];

MaxTimes1 = 10;
interationtimes = 1:10;
a_100_20_100_fix = [a_100_20_100 zeros(1,MaxTimes1-length(a_100_20_100))]; 
a_100_40_100_fix = [a_100_40_100 zeros(1,MaxTimes1-length(a_100_40_100))]; 
a_100_60_100_fix = [a_100_60_100 zeros(1,MaxTimes1-length(a_100_60_100))]; 
a_100_80_100_fix = [a_100_80_100 zeros(1,MaxTimes1-length(a_100_80_100))]; 
a_100_20_300_fix = [a_100_20_300 zeros(1,MaxTimes1-length(a_100_20_300))]; 
a_100_40_300_fix = [a_100_40_300 zeros(1,MaxTimes1-length(a_100_40_300))]; 
a_100_60_300_fix = [a_100_60_300 zeros(1,MaxTimes1-length(a_100_60_300))];  
a_100_80_300_fix = [a_100_80_300 zeros(1,MaxTimes1-length(a_100_80_300))]; 
a_100_20_500_fix = [a_100_20_500 zeros(1,MaxTimes1-length(a_100_20_500))]; 
a_100_40_500_fix = [a_100_40_500 zeros(1,MaxTimes1-length(a_100_40_500))];
a_100_60_500_fix = [a_100_60_500 zeros(1,MaxTimes1-length(a_100_60_500))];
a_100_80_500_fix = [a_100_80_500 zeros(1,MaxTimes1-length(a_100_80_500))];
figure(2)
plot(interationtimes,a_50_10_100_fix,'g-o');  
hold on;  xlabel('迭代的次数'), ylabel('需要重播的比例');
plot(interationtimes,a_50_20_100_fix,'b--o');
plot(interationtimes,a_50_30_100_fix,'m-.o');
plot(interationtimes,a_50_10_300_fix,'g-*');
plot(interationtimes,a_50_20_300_fix,'b--*');
plot(interationtimes,a_50_30_300_fix,'m-.*');
plot(interationtimes,a_50_10_500_fix,'g-+');  
plot(interationtimes,a_50_20_500_fix,'b--+');
plot(interationtimes,a_50_30_500_fix,'m-.+');
legend('QD2D=100,QCellular=10','QD2D=100,QCellular=20','QD2D=100,QCellular=30','QD2D=300,QCellular=10','QD2D=300,QCellular=20','QD2D=300,QCellular=30','QD2D=500,QCellular=10','QD2D=500,QCellular=20','QD2D=500,QCellular=30');


a_150_30_100_fix = [a_150_30_100 zeros(1,MaxTimes1-length(a_150_30_100))]; 
a_150_60_100_fix = [a_150_60_100 zeros(1,MaxTimes1-length(a_150_60_100))]; 
a_150_90_100_fix = [a_150_90_100 zeros(1,MaxTimes1-length(a_150_90_100))]; 
a_150_120_100_fix = [a_150_120_100 zeros(1,MaxTimes1-length(a_150_120_100))]; 
a_150_30_300_fix = [a_150_30_300 zeros(1,MaxTimes1-length(a_150_30_300))]; 
a_150_60_300_fix = [a_150_60_300 zeros(1,MaxTimes1-length(a_150_60_300))]; 
a_150_90_300_fix = [a_150_90_300 zeros(1,MaxTimes1-length(a_150_90_300))];  
a_150_120_300_fix = [a_150_120_300 zeros(1,MaxTimes1-length(a_150_120_300))]; 
a_150_30_500_fix = [a_150_30_500 zeros(1,MaxTimes1-length(a_150_30_500))]; 
a_150_60_500_fix = [a_150_60_500 zeros(1,MaxTimes1-length(a_150_60_500))];
a_150_90_500_fix = [a_150_90_500 zeros(1,MaxTimes1-length(a_150_90_500))];
a_150_120_500_fix = [a_150_120_500 zeros(1,MaxTimes1-length(a_150_120_500))];
figure(3)
plot(interationtimes,a_150_60_100_fix,'g-o');  
hold on;  xlabel('迭代的次数'), ylabel('需要重播的比例');
plot(interationtimes,a_150_90_100_fix,'b--o');
plot(interationtimes,a_150_120_100_fix,'m-.o');
plot(interationtimes,a_150_60_300_fix,'g-*');
plot(interationtimes,a_150_90_300_fix,'b--*');
plot(interationtimes,a_150_120_300_fix,'m-.*');
plot(interationtimes,a_150_60_500_fix,'g-+');  
plot(interationtimes,a_150_90_500_fix,'b--+');
plot(interationtimes,a_150_120_500_fix,'m-.+');
legend('QD2D=100,QCellular=60','QD2D=100,QCellular=90','QD2D=100,QCellular=120','QD2D=300,QCellular=60','QD2D=300,QCellular=90','QD2D=300,QCellular=120','QD2D=500,QCellular=60','QD2D=500,QCellular=90','QD2D=500,QCellular=120');


a_50_10_100_fix = [a_50_10_100 zeros(1,MaxTimes1-length(a_50_10_100))]; 
a_50_20_100_fix = [a_50_20_100 zeros(1,MaxTimes1-length(a_50_20_100))]; 
a_50_30_100_fix = [a_50_30_100 zeros(1,MaxTimes1-length(a_50_30_100))]; 
a_50_40_100_fix = [a_50_40_100 zeros(1,MaxTimes1-length(a_50_40_100))]; 

a_50_10_300_fix = [a_50_10_300 zeros(1,MaxTimes1-length(a_50_10_300))]; 
a_50_20_300_fix = [a_50_20_300 zeros(1,MaxTimes1-length(a_50_20_300))]; 
a_50_30_300_fix = [a_50_30_300 zeros(1,MaxTimes1-length(a_50_30_300))];  
a_50_40_300_fix = [a_50_40_300 zeros(1,MaxTimes1-length(a_50_40_300))]; 

a_50_10_500_fix = [a_50_10_500 zeros(1,MaxTimes1-length(a_50_10_500))]; 
a_50_20_500_fix = [a_50_20_500 zeros(1,MaxTimes1-length(a_50_20_500))];
a_50_30_500_fix = [a_50_30_500 zeros(1,MaxTimes1-length(a_50_30_500))];
a_50_40_500_fix = [a_50_40_500 zeros(1,MaxTimes1-length(a_50_40_500))];

plot(interationtimes,a_50_10_100_fix,'-o');  hold on;  xlabel('迭代的次数')，ylabel('需要重播的比例');
plot(interationtimes,a_50_20_100_fix,':.o');
plot(interationtimes,a_50_30_100_fix,':o');
plot(interationtimes,a_50_10_300_fix,'-*');
plot(interationtimes,a_50_20_300_fix,':.*');
plot(interationtimes,a_50_30_300_fix,':*');
plot(interationtimes,a_50_10_500_fix,'-+');  
plot(interationtimes,a_50_20_500_fix,':.+');
plot(interationtimes,a_50_30_500_fix,':+');


plot(interationtimes,a_50_40_100_fix,'-');  hold on;  xlabel('迭代的次数')，ylabel('需要重播的比例');
plot(interationtimes,a_50_40_300_fix,'-');  hold on;  xlabel('迭代的次数')，ylabel('需要重播的比例');
plot(interationtimes,a_50_40_500_fix,'-');  hold on;  xlabel('迭代的次数')，ylabel('需要重播的比例');
%%%%-------------------------%%%%%
a_100_20_100_fix = [a_100_20_100 
a_100_40_100_fix = [a_100_40_100 
a_100_60_100_fix = [a_100_60_100 
a_100_80_100_fix = [a_100_80_100 

a_100_20_300_fix = [a_100_20_300 
a_100_40_300_fix = [a_100_40_300 
a_100_60_300_fix = [a_100_60_300 
a_100_80_300_fix = [a_100_80_300 

a_100_20_500_fix = [a_100_20_500 
a_100_40_500_fix = [a_100_40_500 
a_100_60_500_fix = [a_100_60_500 
a_100_80_500_fix = [a_100_80_500 

%%%%-------------------------%%%%%
a_150_30_100_fix = [a_150_30_100 
a_150_60_100_fix = [a_150_60_100 
a_150_90_100_fix = [a_150_90_100 
a_150_120_100_fix = [a_150_120_100
a_150_30_300_fix = [a_150_30_300 
a_150_60_300_fix = [a_150_60_300 
a_150_90_300_fix = [a_150_90_300 
a_150_120_300_fix = [a_150_120_300
a_150_30_500_fix = [a_150_30_500 
a_150_60_500_fix = [a_150_60_500 
a_150_90_500_fix = [a_150_90_500 
a_150_120_500_fix = [a_150_120_500 

