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

            % %类型1：cellular之间
            % %所有的蜂窝链接都是相邻的
            % for i = 1:QCellular
            %     for j = 1:QCellular
            %         if i~=j
            %             Vtrue{i} = [Vtrue{i} j];
            %         end
            %     end
            % end

            %类型2：D2D受到cellular干扰
            for i = 1:QCellular
                for j = QCellular+1:Q
                    if (cellular_loc(i,1)-X_Yr(j-QCellular,1))^2+(cellular_loc(i,2)-X_Yr(j-QCellular,2))^2<(d_broadcast(i))^2 &&i~=j
                        Vtrue{i} = [Vtrue{i} j];
                    end
                end
            end

            %类型3：celllular受到D2D之间

            for i = QCellular+1:Q          
                if D2D_loc(i-QCellular,1)+D2D_loc(i-QCellular,2)^2 < (d_broadcast(i))^2 && i~=j
                    Vtrue{i} = [Vtrue{i} 10000]; %用一个特殊的数字表示基站
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%     stage1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            flag = ones(1,Q);    %标记变量，初始化全1
            Nq = 0;              %统计迭代次数
            storage = [];        %存放每次迭代pp的值，行向量

            %-------   初始化V_nbor(每个D2D链路广播范围内的其他D2D链路)   -------%
            
            V_nbor_all = Vtrue; %V_nbor_all 用于存放每个D2D链路能收到广播范围内的【所有】D2D和cellular link
            V_nbor = cell(Q,1); %V_nbor 用于存放每个链路能收到广播范围内的D2Dlink和cellular（不断变化）
            for i =1:Q          %开始的时候全部能收到
                V_nbor{i} = V_nbor_all{i};
            end

            while (sum(flag)~=0)
        
                link = find(flag==1);   %flag = 1的link,行向量
                
                %---------------------      分配资源    -----------------------%
                RU = zeros(1,Q);   %行向量,每条flag=1链路的资源分配情况
                for i = 1:length(link)
                    RU(1,link(1,i)) = randint(1,1,[1,N]);
                end
                
                link0 = find(RU == 0);      %未分配资源的link
                [Ulink,U1] = cw_stage1(RU,link0); %Ulink为重复的RU的link，U1为重复RU（去0）

                K = setdiff(RU,[U1 0]);     %不重复的RU，为行向量
                                            %setdiff函数为求两个元素的差，相当于 集合K = 集合RU - 集合[U1 0]

                [tf H] = ismember(K,RU);    %H返回不重复RU的link，为行向量，与K同维
                                            %tf返回一个0、1集合，如果元素在RU里面则为1，否则为0  


                %--------------  求重复的link和RU（邻居范围内的）  -------------%
                V_nbor_RU = cell(Q,1);      %广播范围内的link所分配的RU
                for i = 1:Q
                    V_nbor_RU{i} = [];
                    for j = 1:length(V_nbor{i})
                        if V_nbor{i}(j)==10000
                            V_nbor_RU{i} = [V_nbor_RU{i} RU(link(1:QCellular))];
                        else   
                            V_nbor_RU{i} = [V_nbor_RU{i} RU(V_nbor{i}(j))];
                        end
                        
                    end
                end
                
                U_nbor = cell(Q,1);    %重复的RU,即Uq(供stage2使用)
                for i = 1:Q
                    if ~isempty(V_nbor_RU{i})
                        b = tabulate(V_nbor_RU{i});
                        c = b(find(b(:,2)>1),1);
                        U_nbor{i} = setdiff(c',0);     %行向量
                    else
                        U_nbor{i} = [];
                    end
                end
                
                %----------------        初步构建干扰图       -----------------%
                for i = 1:Q
                    V{i} = intersect(V_nbor{i},H);     %查用法是否无误
                end
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%     stage2     %%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %----------------------      分配时延    ----------------------%
                AckDelay = [];
                AckDelay = unidrnd(M,1,Q);  %给每个D2D链路均分配一个随机回退的延时用于发送ACK协议
                
                [mk] = cw_stage2(AckDelay);   %mu为被分配多个link的时延,mu_link为对应的link    
                [tf1 mk_link] = ismember(mk,AckDelay);  %mk_link为只分配一个link的时延对应的link
                
                Ack_nbor_Delay = cell(Q,1);   %邻居分配的时延
                for i = 1:Q
                    Ack_nbor_Delay{i} = [];   %初始化为空
                    for j = 1:length(V_nbor{i})
                        if V_nbor{i}(j)==10000
                            Ack_nbor_Delay{i} = [Ack_nbor_Delay{i} AckDelay(link(1:QCellular))];
                        else
                            Ack_nbor_Delay{i} = [Ack_nbor_Delay{i} AckDelay(V_nbor{i}(j))];
                        end
                    end
                end
                
                V_Delay_Fenlei = cell(Q,1);   %将一个延时有多个link的归类（cell套cell）
                for i =1:Q
                    if ~isempty(Ack_nbor_Delay{i})
                        b = tabulate(Ack_nbor_Delay{i});
                        c = b(find(b(:,2)>1),1);
                        c = c';     %行向量
                        ca = cell(length(c),1);
                        for j = 1:length(c)
                            d = [];
                            d = intersect(find(AckDelay == c(j)),V_nbor{i});
                            ca{j} = [c(j) d];
                        end
                        V_Delay_Fenlei{i} = ca;
                    else
                        V_Delay_Fenlei{i} = [];
                    end
                end
                
                %------------------     相同时延下的ACK    --------------------%
                V_Delay_Fenlei_U = cell(Q,1);       %同一时延下发送的U的集合（存放RU）（cell套cell）
                V_Delay_Fenlei_UU = cell(Q,1);      %同一时延下发送的重复的U的集合（cell套cell）
                V_Delay_Fenlei_Ulink = cell(Q,1);   %同一时延下发送的重复的U的集合的link（cell套cell）
                V_Delay_Fenlei_uu = cell(Q,1);      %所有重复U的集合
                V_Delay_Fenlei_ulink = cell(Q,1);   %所有重复U的集合的link
                
                for i = 1:Q      %V_Delay_Fenlei_U
                    for j = 1:size(V_Delay_Fenlei{i},1)
                        V_Delay_Fenlei_U{i}{j,1} = [];
                        for k = 2:length(V_Delay_Fenlei{i}{j,1})
                            V_Delay_Fenlei_U{i}{j,1} = [V_Delay_Fenlei_U{i}{j,1} U_nbor{(V_Delay_Fenlei{i}{j,1}(k))}];
                        end
                    end     %到此循环，V_Delay_Fenlei_U{i}就得到了
                end
                
                for i = 1:Q   %V_Delay_Fenlei_UU
                    for j = 1:size(V_Delay_Fenlei_U{i},1)
                        if ~isempty(V_Delay_Fenlei_U{i}{j})
                            ua = tabulate(V_Delay_Fenlei_U{i}{j,1});
                            ub = ua(find(ua(:,2)>1),1);
                            V_Delay_Fenlei_UU{i}{j,1} = ub';   %行向量
                        end
                    end
                end
                
                for i = 1:Q    %V_Delay_Fenlei_Ulink
                    for j = 1:size(V_Delay_Fenlei_UU{i},1)
                        ux = [];     %临时存储link
                        for k = 1:length(V_Delay_Fenlei_UU{i}{j,1})
                            ux = [ux find(RU == V_Delay_Fenlei_UU{i}{j,1}(k))];
                        end
                        V_Delay_Fenlei_Ulink{i}{j,1} = intersect(V_nbor{i},ux);   %可能邻居之外的link也使用同样的RU，故求交集
                    end
                end
                
                for i = 1:Q   %V_Delay_Fenlei_uu
                    V_Delay_Fenlei_uu{i}= [];
                    for j = 1:size(V_Delay_Fenlei_UU{i},1)
                        V_Delay_Fenlei_uu{i} = [V_Delay_Fenlei_uu{i} V_Delay_Fenlei_UU{i}{j}];  %将UU中的cell变成一个行向量
                    end
                    V_Delay_Fenlei_uu{i} = unique(V_Delay_Fenlei_uu{i});
                end
                
                for i = 1:Q   %V_Delay_Fenlei_ulink
                    V_Delay_Fenlei_ulink{i} = [];
                    for j = 1:size(V_Delay_Fenlei_Ulink{i},1)
                        V_Delay_Fenlei_ulink{i} = [V_Delay_Fenlei_ulink{i} V_Delay_Fenlei_Ulink{i}{j}];  %将UU中的cell变成一个行向量
                    end
                    V_Delay_Fenlei_ulink{i} = unique(V_Delay_Fenlei_ulink{i});
                end
                
                
                %-----------------------    更新干扰图   ----------------------%
                for i = 1:Q
                    V{i} = setdiff(V_nbor_all{i},V_Delay_Fenlei_ulink{i});
                end
                %----------------------    更新flag     ----------------------%
                count = zeros(1,Q);
                for i = 1:Q
                    if flag(i) ~=0 && length(setdiff(V_Delay_Fenlei_uu{i},RU(i)))~= length(V_Delay_Fenlei_uu{i})  %RU(i)在V_Delay_Fenlei_ulink{i}中
                        count(i) = 1;
                    end
                end
                flag = count;
                
                %----------------------  更新V_nbor    -----------------------%
                count_link = find(count==1);    %flag=1(即count=1)的link
                for i =1:Q
                    V_nbor{i} = intersect(V_nbor{i},count_link);
                end
                
                Nq = Nq+1;
                pp = length(count_link)/Q;    %flag=1的link个数/总D2Dlink数
                storage = [storage pp];
            

            end
            data{x,y} = storage;
              

        end
    end
end
toc;
