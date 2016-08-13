%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%           试验8     算法1（协议）    M=2     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
tic;
%load data_of_link;                   %'D2D_loc','d_broadcast','X_Yr','D','b','arr_frame','lev_frame'
r_cell = 300;                        %小区半径
QD2D = 500;                          %D2D链路对数
QCellular = 30;                      %cellular链路对数
N = 50;                              %资源块个数
M = 2;                               %最大回退延时
Dmax = 50;                           %D2Dlink最大通信半径
u_frame_D2D = 30;                         %每frame进入/离开的D2D用户个数
u_frame_cellular = 3;                      %每frame进入/离开的cellular link用户个数
w = 2;
maxframe = 400;                     %总帧

Q = QD2D + QCellular;
V = cell(Q,1);                       %干扰图
Vtrue = cell(Q,1);                  %真实干扰图1帧
data_pp_M2 = [];                        %存储结果

%%%%%%%%----------          frame=0时候cellular link位置           ----------%%%%%%%%
cellular_loc = unifrnd(-r_cell,r_cell,QCellular,2);     %小区内每个D2D链路的位置
for i = 1:QCellular
    while sum(cellular_loc(i,:).^2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
        cellular_loc(i,1)= rand(1)*2*r_cell-r_cell;
        cellular_loc(i,2)= rand(1)*2*r_cell-r_cell;
    end
end
Dcellular = sqrt(sum(cellular_loc.^2,2));      %发射端与接收端距离

%%%%%%%%----------          frame=0时候D2D link位置           ----------%%%%%%%%
D2D_loc = unifrnd(-r_cell,r_cell,QD2D,2);    %小区内每个D2D链路的位置
for i = 1:QD2D
    while sum(D2D_loc(i,:).^2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
        D2D_loc(i,1)= rand(1)*2*r_cell-r_cell;
        D2D_loc(i,2)= rand(1)*2*r_cell-r_cell;
    end
end
D_D2D = unifrnd(1,Dmax,QD2D,1);      %发射端与接收端距离
A_D2D = unifrnd(0,2*pi,QD2D,1);      %极坐标角度在0到2pi之间

D = [Dcellular',D_D2D'];
d_broadcast = D * w;

X_Yr=zeros(QD2D,2);     %接收端坐标
for i = 1:QD2D
    X_Yr(i,1)=D2D_loc(i,1)+D(i)*cos( A_D2D(i));
    X_Yr(i,2)=D2D_loc(i,2)+D(i)*sin( A_D2D(i));
end

%%%%%%---------          新进入cellular link用户的进入帧，位置            ----------%%%%%%

arr_loc_cellular = unifrnd(-r_cell,r_cell,(u_frame_cellular+1)*(maxframe+1),2); %新用户位置
for i = 1:(u_frame_cellular+1)*(maxframe+1)
    while sum(arr_loc_cellular(i,:).^2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
        arr_loc_cellular(i,1)= rand(1)*2*r_cell-r_cell;
        arr_loc_cellular(i,2)= rand(1)*2*r_cell-r_cell;
    end
end

%%%%%%---------          新进入D2D用户的进入帧，位置            ----------%%%%%%
arr_loc_D2D = unifrnd(-r_cell,r_cell,(u_frame_D2D+1)*(maxframe+1),2); %新用户位置
for i = 1:(u_frame_D2D+1)*(maxframe+1)
    while sum(arr_loc_D2D(i,:).^2) > r_cell*r_cell   %如果随机生成的发射端在圆之外，就继续生成，直到在圆之内为止
        arr_loc_D2D(i,1)= rand(1)*2*r_cell-r_cell;
        arr_loc_D2D(i,2)= rand(1)*2*r_cell-r_cell;
    end
end

%%%--------                      变量初始化                     ---------%%%
d_broadcast_1frame = d_broadcast;
X_Yr_1frame = X_Yr;
cellular_loc_1frame = cellular_loc;
D2D_loc_1frame = D2D_loc;
QD2D_temp = QD2D;
QCellular_temp = QCellular;
Q_temp = QD2D_temp + QCellular_temp;

% link_types = [ones(1,QCellular),zeros(1,QD2D)];  %标记链接的类型
% link_types_1frame = link_types;

d_broadcast_cellular = d_broadcast(1:QCellular);
d_broadcast_D2D = d_broadcast(QCellular+1:Q);
d_broadcast_1frame_cellular = d_broadcast_cellular;
d_broadcast_1frame_D2D = d_broadcast_D2D;





flag = ones(1,Q);
count = zeros(1,Q);

for slot = 0:maxframe-1
    %%----------------          每1帧新进来的cellular link         ----------------%%
    loc_new_cellular = zeros(u_frame_cellular,2);
    for i = 1:u_frame_cellular
        loc_new_cellular(i,:) = arr_loc_cellular((slot+1)*u_frame_cellular+i,:);    %n*2的矩阵，坐标为两个值
    end
    loc_temp_cellular = [D2D_loc_1frame;loc_new_cellular];
    D_arr_cellular = sqrt(sum(loc_new_cellular.^2,2));    %cellular link 与基站的距离
    Bro_add_arr_cellular = D_arr_cellular * w;                             %增加link的广播半径
    d_broadcast_1frame_cellular = [d_broadcast_1frame_cellular,Bro_add_arr_cellular'];
    % link_types = [link_types,ones(1,length(Bro_add_arr_cellular'))]; %标记增加的类型

    flag = [flag(1:QCellular) ones(1,u_frame_cellular) flag(QCellular+1:Q)];
    QCellular = QCellular_temp + u_frame_cellular;
    % flag = [flag ones(1,u_frame_cellular)];

    %%----------------          每1帧新进来的D2D link         ----------------%%
    loc_new_D2D = zeros(u_frame_D2D,2);
    for i = 1:u_frame_D2D
        loc_new_D2D(i,:) = arr_loc_D2D((slot+1)*u_frame_D2D+i,:);    %n*2的矩阵，坐标为两个值
    end
    loc_temp_D2D = [D2D_loc_1frame;loc_new_D2D];

    D_arr_D2D = unifrnd(1,Dmax,u_frame_D2D,1);      %发射端与接收端距离
    A_arr_D2D = unifrnd(0,2*pi,u_frame_D2D,1);      %极坐标角度在0到2pi之间
    Bro_add_arr_D2D = D_arr_D2D * w;                             %增加link的广播半径
    d_broadcast_1frame_D2D = [d_broadcast_1frame_D2D,Bro_add_arr_D2D'];
    % link_types = [link_types,zeros(1,length(Bro_add_arr_D2D))]; %标记增加的类型

    XYr_arr_1frame = zeros(u_frame_D2D,2);
    for i = 1:u_frame_D2D
        XYr_arr_1frame(i,1) = loc_new_D2D(i,1)+D_arr_D2D(i)*cos(A_arr_D2D(i));
        XYr_arr_1frame(i,2) = loc_new_D2D(i,2)+D_arr_D2D(i)*sin(A_arr_D2D(i));
    end
    XYr_temp = [X_Yr_1frame;XYr_arr_1frame];    %接收端坐标
    QD2D = QD2D_temp + u_frame_D2D;
    flag = [flag ones(1,u_frame_D2D)];

    Q = QD2D + QCellular;

    %%-----------             每1帧离开的cellular link             ----------------%%
    lev_num_cellular = u_frame_cellular;
    all_link_cellular = randperm(QCellular);
    lev_link_cellular_unsort = all_link_cellular(1:u_frame_cellular);
    lev_link_cellular = sort(lev_link_cellular_unsort);

    lev_link_cellular_index =[];
    for i = 1:u_frame_cellular
        lev_link_cellular_index =[lev_link_cellular_index,lev_link_cellular(u_frame_cellular-i+1)];
    end
    loc_temp_cellular(lev_link_cellular_index,:) = [];
    d_broadcast_1frame_cellular(lev_link_cellular_index) = [];
    % flag(lev_link_cellular_index) =[];

    %%-----------             每1帧离开的D2D link             ----------------%%
    lev_num_D2D = u_frame_D2D;
    all_link_D2D = randperm(QD2D);
    lev_link_D2D_unsort = all_link_D2D(1:u_frame_D2D);
    lev_link_D2D = sort(lev_link_D2D_unsort);

    lev_link_D2D_index =[];
    for i = 1:u_frame_D2D
        lev_link_D2D_index =[lev_link_D2D_index,lev_link_D2D(u_frame_D2D-i+1)];
    end
    loc_temp_D2D(lev_link_D2D_index,:) = [];
    XYr_temp(lev_link_D2D_index,:) = [];
    d_broadcast_1frame_D2D(lev_link_D2D_index) = [];

    lev_link_flag_index =[lev_link_cellular_index,lev_link_D2D_index + QCellular];
    flag(lev_link_flag_index) =[];


    D2D_loc_1frame = loc_temp_D2D;
    cellular_loc_1frame = loc_temp_cellular;
    X_Yr_1frame = XYr_temp;

    d_broadcast_1frame = [d_broadcast_1frame_cellular,d_broadcast_1frame_D2D];

    QCellular = QCellular - lev_num_cellular;
    QD2D = QD2D - lev_num_D2D;
    Q = QD2D + QCellular;

    %-------------------          构建真实干扰图        -------------------%
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
    %--------------     统计边数      --------------%
    num_true = 0;
    for i = 1 : Q
        num_true = num_true + length(Vtrue{i});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%     stage1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    link_tpyes = [ones(1,QCellular),zeros(1,QD2D)];  %标记链接的类型
    cellular_link_type = find(link_tpyes==1);
    D2D_link_tpye = find(link_tpyes==0);

    %-------   初始化V_nbor(每个D2D链路广播范围内的其他D2D链路)   -------%
    V_nbor_all = Vtrue; %V_nbor用于存放每个D2D链路能收到广播范围内的【所有】D2Dlink
    V_nbor = cell(Q,1); %V_nbor用于存放每个D2D链路能收到广播范围内的D2Dlink（不断变化）
    for i =1:Q
        V_nbor{i} = V_nbor_all{i};
    end
    link = find(flag==1);   %flag = 1的link,行向量

    cellular_link = intersect(cellular_link_type,link);

    %---------------------      分配资源    -----------------------%
    RU = zeros(1,Q);   %行向量,每条flag=1链路的资源分配情况
    for i = 1:length(link)
        RU(1,link(1,i)) = randi([1,N],1,1);
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
                V_nbor_RU{i} = [V_nbor_RU{i} RU(cellular_link)];
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
    V1 = cell(Q,1);
    for i = 1:Q
        V1{i} = intersect(V_nbor{i},H);     %查用法是否无误
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
                Ack_nbor_Delay{i} = [Ack_nbor_Delay{i} AckDelay(cellular_link)];
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
        V{i} = union(setdiff(V_nbor_all{i},V_Delay_Fenlei_ulink{i}),V1{i});
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
    num_simu = 0;    %协议干扰边的个数
    for i = 1:Q
        num_simu = num_simu +length(V{i});
    end
    Q_temp1 = Q;
    mis_edge = abs(num_true-num_simu);
    pp = mis_edge/num_true;
    data_pp_M2 = [data_pp_M2 pp]
    slot = slot + 1;
    save data_pp_M2_30,'data_pp_M2';

end

toc;
