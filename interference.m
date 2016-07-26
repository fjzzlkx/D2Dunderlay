%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        干扰图构建协议   M=5      %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
tic;
r_cell = 300;                        %小区半径
QD2Dpool = [100 300 500];            %D2D链路对数

QCellularpool = [10 20 30 40;
                 20 40 60 80;
                 30 60 90 120];        %蜂窝链接的对

Npool = [50 100 150];                %资源块个数
M = 5;                               %最大回退延时
K0 = 3;                              %目标信干比
S0 = 0;                              %定义广播半径的信干比阈值
sigma = 8;                           %衰落标准差
Dmax = 50;                           %D2Dlink最大通信半径


data = cell(size(QCellularpool,1)*size(QCellularpool,2), length(Npool));   %存放storage
for x = 1:length(Npool)
    for y = 1:length(QD2Dpool)
        for z = 1:size(QCellularpool,2)
            QD2D = QD2Dpool(y);
            QCellular = QCellularpool(y,z);
            Q = QD2D + QCellular;
            N = Npool(x);
            V = cell(Q,1);      %干扰图
            Vtrue = cell(Q,1);  %真实干扰图

            %%%%---------       cellular link在整个小区内均匀分布       ---------%%%%%
            cellular_loc = unifrnd(-r_cell,r_cell,QCellular,2);    %小区内每个D2D链路的位置
            for i = 1:QCellular
                while cellular_loc(i,1)*cellular_loc(i,1)+cellular_loc(i,2)*cellular_loc(i,2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
                    cellular_loc(i,1)= rand(1)*2*r_cell-r_cell;
                    cellular_loc(i,2)= rand(1)*2*r_cell-r_cell;
                end
                % plot( cellular_loc(i,1), cellular_loc(i,2),'+')
                % hold on;
            end
            Dcellular = sqrt(sum(cellular_loc.^2,2));      %发射端与接收端距离
            %Acellular = unifrnd(0,2*pi,QCellular,1);      %极坐标角度在0到2pi之间



            d_broadcast = zeros(Q,1);   %存放所有的广播半径


            for i = 1:QCellular         %存放蜂窝的广播半径
                d_broadcast(i) = Dcellular(i)*2;
            end
            
            %%%%-----------------D2D连接实在圆形区域内均匀分布-------------------%%%%
            D2D_loc = unifrnd(-r_cell,r_cell,QD2D,2);    %小区内每个D2D链路的位置
            for i = 1:QD2D
                while D2D_loc(i,1)*D2D_loc(i,1)+D2D_loc(i,2)*D2D_loc(i,2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
                    D2D_loc(i,1)= rand(1)*2*r_cell-r_cell;
                    D2D_loc(i,2)= rand(1)*2*r_cell-  r_cell;
                end
            end
            DD2D = unifrnd(1,Dmax,QD2D,1);      %发射端与接收端距离
            AD2D = unifrnd(0,2*pi,QD2D,1);         %极坐标角度在0到2pi之间

            for i = 1:QD2D                      %存放D2D链接的广播半径
                d_broadcast(i+QCellular) = DD2D(i)*2;
            end
            
            X_Yr=zeros(QD2D,2);     %接收端坐标
            for i = 1:QD2D         
                X_Yr(i,1)=D2D_loc(i,1)+DD2D(i)*cos( AD2D(i));
                X_Yr(i,2)=D2D_loc(i,2)+DD2D(i)*sin( AD2D(i));
            end

            %-------------------       构建真实干扰图       -------------------%   
            for i = 1:Q      %初始化Vtrue
                Vtrue{i} = [];
            end

            %类型1：cellular之间
            %所有的蜂窝链接都是相邻的
            for i = 1:QCellular
                for j = 1:QCellular
                    if i~=j
                        Vtrue{i} = [Vtrue{i} j];
                    end
                end
            end

            %类型2：D2D受到cellular干扰
            for i = 1:QCellular
                for j = QCellular+1:Q
                    if (cellular_loc(i,1)-X_Yr(j-QCellular,1))^2+(cellular_loc(i,2)-X_Yr(j-QCellular,2))^2<(d_broadcast(i))^2 &&i~=j
                        Vtrue{i} = [Vtrue{i} j];
                    end
                end
            end

            %类型3：celllular受到D2D干扰

            for i = QCellular+1:Q          
                if D2D_loc(i-QCellular,1)+D2D_loc(i-QCellular,2)^2 < (d_broadcast(i))^2 && i~=j
                    %for j = 1:QCellular
                    %   Vtrue{i} = [Vtrue{i} j];
                    %end
                    Vtrue{i} = [Vtrue{i} 10000];%选取一个特殊的数字当作基站
                end
            end

            %类型4：D2D之间
            for i = QCellular+1:Q
                for j = QCellular+1:Q
                    if (D2D_loc(i-QCellular,1)-X_Yr(j-QCellular,1))^2+(D2D_loc(i-QCellular,2)-X_Yr(j-QCellular,2))^2<(d_broadcast(i))^2 &&i~=j
                        Vtrue{i} = [Vtrue{i} j];
                    end
                end
            end
            

        end
    end
end
toc;
