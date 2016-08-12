%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%           试验8     算法1（协议）    M=2     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
tic;
load data_of_link;                   %'D2D_loc','d_broadcast','X_Yr','D','b','arr_frame','lev_frame'
r_cell = 300;                        %小区半径
Q_D2D = 700;                         %D2D链路对数
Q_cellular = 30;;                    %cellular链路对数                     
N = 50;                              %资源块个数
M = 2;                               %最大回退延时
Dmax = 50;                           %D2Dlink最大通信半径
u_D2D = 200;                         %每秒进入/离开的D2D用户个数
u_cellular = 20;                     %每秒进入/离开的cellular link用户个数
w = 2;
maxframe = 1200;                     %总帧
QCellular_arr = maxframe*u_cellular / 100;              %新加入cellular link用户的总个数
QCellular_lev = QCellular_arr;                          %离开cellular link用户的总个数
QD2D_arr = maxframe* u_D2D / 100;               %新加入D2D用户的总个数
QD2D_lev = QD2D_arr;                            %离开D2D用户的总个数



V = cell(Q,1);                       %干扰图
Vtrue = cell(Q,1);                   %真实干扰图(M+1)帧
Vtrue1 = cell(Q,1);                  %真实干扰图1帧
data_pp = [];                        %存储结果

%%%--------                      变量初始化                     ---------%%%
d_broadcast1 = d_broadcast;
X_Yr1 = X_Yr;
cellular_loc1 = cellular_loc;
D2D_loc1 = D2D_loc;
slot = 0;
QD2D_temp = QD2D;
QD2D_temp1 = QD2D;
QCellular_temp = QCellular;
QCellular_temp1 = QCellular;
Q_temp = QD2D_temp + QCellular_temp;
Q_temp1 = QD2D_temp1 + QCellular_temp1;


link_tpyes = [ones(1,QCellular),zeros(1,QD2D)];  %标记链接的类型
% cellular_link_type = find(link_tpyes==1);
% D2D_link_tpye = find(link_tpyes==0);


while (slot <= maxframe-M-1)
    %%----------------       每(M+1)帧新进来的cellular link       ----------------%%
    arrx_cellular = find(arr_frame_cellular >= (slot));
    arry_cellular = find(arr_frame_cellular <= slot+M+1);
    arr_cellular_temp = intersect(arrx_cellular,arry_cellular);  %当前新加入的link的位置
    loc_new_cellular = zeros(length(arr_cellular_temp),2);    %n*2的矩阵，坐标为两个值
    
    if length(arr_cellular_temp) ~= 0
        for i = 1:length(arr_cellular_temp)
            loc_new_cellular(i,1) = arr_loc_cellular(arr_cellular_temp(i),1);
            loc_new_cellular(i,2) = arr_loc_cellular(arr_cellular_temp(i),2);
        end
        QCellular_temp = QCellular_temp + length(arr_cellular_temp);
        loc_cellular_temp = [cellular_loc;loc_new_cellular];

        %D_arr = unifrnd(1,Dmax,length(arr_temp),1);       
        D_arr_cellular = sqrt(sum(arr_cellular_temp.^2,2));    %cellular link 与基站的距离
        
        Bro_add_arr_cellular = D_arr_cellular * w;                             %增加link的广播半径
        d_broadcast = [d_broadcast;Bro_add_arr_cellular];
        
        link_tpyes = [link_tpyes,ones(length(Bro_add_arr_cellular))]; %标记增加的类型
    else
        loc_cellular_temp = cellular_loc;
        d_broadcast = d_broadcast;
        link_tpyes =link_tpyes;
    end
    
    %%----------------       每(M+1)帧新进来的D2Dlink        ----------------%%
    arrx_D2D = find(arr_frame_D2D >= (slot));
    arry_D2D = find(arr_frame_D2D <= slot+M+1);
    arr_D2D_temp = intersect(arrx_D2D,arry_D2D);  %当前新加入的link的位置
    loc_new_D2D = zeros(length(arr_D2D_temp),2);    %n*2的矩阵，坐标为两个值
    if length(arr_D2D_temp) ~= 0
        for i = 1:length(arr_D2D_temp)
            loc_new_D2D(i,1) = arr_loc_D2D(arr_D2D_temp(i),1);
            loc_new_D2D(i,2) = arr_loc_D2D(arr_D2D_temp(i),2);
        end
        QD2D_temp = QD2D_temp + length(arr_D2D_temp);
        loc_D2D_temp = [D2D_loc;loc_new_D2D];

        XYr_arr=zeros(length(arr_D2D_temp),2);                   %新增加link的接收端坐标
        D_arr_D2D = unifrnd(1,Dmax,length(arr_D2D_temp),1);      %发射端与接收端距离
        a_arr_D2D = unifrnd(0,2*pi,length(arr_D2D_temp),1);      %极坐标角度在0到2pi之间
        Bro_add_arr_D2D = D_arr_D2D * w;                         %增加D2D link的广播半径
        d_broadcast = [d_broadcast;Bro_add_arr_D2D];
        link_tpyes = [link_tpyes,zeros(length(Bro_add_arr_cellular))]; %标记增加的类型
        
        for i = 1:length(arr_D2D_temp)
            XYr_arr(i,1) = loc_new_D2D(i,1)+D_arr_D2D(i)*cos(a_arr_D2D(i));
            XYr_arr(i,2) = loc_new_D2D(i,2)+D_arr_D2D(i)*sin(a_arr_D2D(i));
        end
        XYr_temp = [X_Yr;XYr_arr];    %接收端坐标
    else
        XYr_temp = X_Yr;
        loc_D2D_temp = D2D_loc;
        d_broadcast = d_broadcast;
        link_tpyes =link_tpyes;
    end

    %%-----------           每(M+1)帧离开的cellular link          ----------------%%
    levx_cellular = find(lev_frame_cellular >= (slot));
    levy_cellular = find(lev_frame_cellular <= slot+M+1);
    lev_num_cellular = length(intersect(levx_cellular,levy_cellular));  %离开cellular link的个数
    if length(lev_num_cellular) ~= 0
        lev_link_cellular = randi([1 QCellular_temp],1,lev_num_cellular);
        for i = 1:length(lev_link_cellular)
            %loc_temp(lev_link_cellular(i),:) = [];
            %XYr_temp(lev_link_cellular(i),:) = [];
            loc_cellular_temp(lev_link_cellular(i),:) = [];
            cellular_link = find(link_tpyes == 1);
            tempIndex1 = cellular_link(lev_link_cellular(i));
            d_broadcast(tempIndex1) = [];
            link_tpyes(tempIndex1) = [];
        end
        cellular_loc = loc_cellular_temp;
        d_broadcast = d_broadcast;
        QCellular_true = QCellular_temp - lev_num_cellular;
        link_tpyes = link_tpyes;
        
    else
        QCellular_true = QCellular_temp;    %求当下真实干扰图需要更新的变量（4个）
        cellular_loc = loc_cellular_temp;
        d_broadcast = d_broadcast;
        link_tpyes = link_tpyes;
    end
    %%-----------           每(M+1)帧离开的D2D link          ----------------%%
    levx_D2D = find(lev_frame_D2D >= (slot));
    levy_D2D = find(lev_frame_D2D <= slot+M+1);
    lev_num_D2D = length(intersect(levx_D2D,levy_D2D));  %离开link的个数
    if length(lev_num_D2D) ~= 0
        lev_link_D2D = randi([1 Q_temp],1,lev_num_D2D);
        for i = 1:length(lev_link_D2D)
            loc_D2D_temp(lev_link_D2D(i),:) = [];
            XYr_temp(lev_link_D2D(i),:) = [];
            D2D_link = find(link_tpyes == 0);
            tempIndex2 = D2D_link(lev_link_D2D(i));
            d_broadcast(tempIndex2) = [];
        end
        D2D_loc = loc_D2D_temp;
        X_Yr = XYr_temp;
        d_broadcast = d_broadcast;
        QD2D_true = QD2D_temp - lev_num_D2D;
        link_tpyes = link_tpyes;
    else
        QD2D_true = QD2D_temp - lev_num_D2D;    %求当下真实干扰图需要更新的变量（5个）
        D2D_loc = loc_D2D_temp;
        X_Yr = XYr_temp;
        d_broadcast = d_broadcast;
        link_tpyes = link_tpyes;
    end


    
    %-------------------         构建真实干扰图         -------------------%
    Q_true = QD2D_true + QCellular_true;
    for i = 1:Q_true      %初始化Vtrue
        Vtrue{i} = [];
    end
    d_broadcast_cellular = d_broadcast(find(d_broadcast == 1));
    d_broadcast_D2D = d_broadcast(find(d_broadcast == 1));
    d_broadcast_sort = [d_broadcast_cellular,d_broadcast_D2D];
    %类型1：cellular之间
    %所有的蜂窝链接都是相邻的
    for i = 1:QCellular_true
        for j = 1:QCellular_true
            if i~=j
                Vtrue{i} = [Vtrue{i} j];
            end
        end
    end

    %类型2：D2D受到cellular干扰
    for i = 1:QCellular_true
        for j = QCellular_true+1:Q
            if (cellular_loc(i,1)-X_Yr(j-QCellular_true,1))^2+(cellular_loc(i,2)-X_Yr(j-QCellular_true,2))^2<(d_broadcast_sort(i))^2 &&i~=j
                Vtrue{i} = [Vtrue{i} j];
            end
        end
    end

    %类型3：celllular受到D2D之间

    for i = QCellular+1:Q
        if D2D_loc(i-QCellular,1)+D2D_loc(i-QCellular,2)^2 < (d_broadcast_sort(i))^2 && i~=j
            Vtrue{i} = [Vtrue{i} 10000]; %用一个特殊的数字表示基站
        end
    end

    %类型4：D2D之间
    for i = QCellular+1:Q
        for j = QCellular+1:Q
            if (D2D_loc(i-QCellular,1)-X_Yr(j-QCellular,1))^2+(D2D_loc(i-QCellular,2)-X_Yr(j-QCellular,2))^2<(d_broadcast_sort(i))^2 &&i~=j
                Vtrue{i} = [Vtrue{i} j];
            end
        end
    end
    num_true = 0;     %真实干扰边的个数
    for i = 1:Q_true
        num_true = num_true + length(Vtrue{i});
    end

    %%----------------          每1帧新进来的link         ----------------%%
    arr_x = find(arr_frame>=(slot));
    arr_y = find(arr_frame<=slot+1);
    atemp = intersect(arr_x,arr_y);  %当前新加入的link的位置
    loc_new_D2D1 = zeros(length(atemp),2);    %n*2的矩阵，坐标为两个值
    if length(atemp) ~= 0
        for i = 1:length(atemp)
            loc_new_D2D1(i,1) = arr_loc(atemp(i),1);
            loc_new_D2D1(i,2) = arr_loc(atemp(i),2);
        end
        Q_temp1 = Q_temp1 + length(atemp);
        loc_temp1 = [D2D_loc1;loc_new_D2D1];

        XYr_arr1=zeros(length(atemp),2);               %新增加link的接收端坐标
        D_arr1 = unifrnd(1,Dmax,length(atemp),1);      %发射端与接收端距离
        b_arr1 = unifrnd(0,2*pi,length(atemp),1);      %极坐标角度在0到2pi之间
        Bro_add1 = D_arr1 * w;                         %增加link的广播半径
        d_broadcast1 = [d_broadcast1;Bro_add1];
        for i = 1:length(atemp)
            XYr_arr1(i,1) = loc_new_D2D1(i,1)+D_arr1(i)*cos( b_arr1(i));
            XYr_arr1(i,2) = loc_new_D2D1(i,2)+D_arr1(i)*sin( b_arr1(i));
        end
        XYr_temp1 = [X_Yr1;XYr_arr1];    %接收端坐标
    else
        XYr_temp1 = X_Yr1;
        loc_temp1 = D2D_loc1;
        d_broadcast1 = d_broadcast1;
    end

    %%-----------             每1帧离开的link             ----------------%%
    levx1 = find(lev_frame>=(slot));
    levy1 = find(lev_frame<=slot+1);
    lev_num1 = length(intersect(levx1,levy1));  %离开link的个数
    if length(lev_num1) ~= 0
        lev_link1 = randint(1,lev_num1,[1 Q_temp1]);
        for i = 1:length(lev_link1)
            loc_temp1(lev_link1(i),:) = [];
            XYr_temp1(lev_link1(i),:) = [];
            d_broadcast1(lev_link1(i)) = [];
        end
        D2D_loc1 = loc_temp1;
        X_Yr1 = XYr_temp1;
        d_broadcast1 = d_broadcast1;
        Q_true1 = Q_temp1 - lev_num1;
    else
        Q_true1 = Q_temp1;    %求当下真实干扰图需要更新的变量（3个）
        D2D_loc1 = loc_temp1;
        X_Yr1 = XYr_temp1;
        d_broadcast1 = d_broadcast1;
    end
    %-------------------          构建真实干扰图        -------------------%
    for i = 1:Q_true1
        Vtrue1{i} = [];
        for j = 1:Q_true1
            if (D2D_loc1(i,1)-X_Yr1(j,1))^2+(D2D_loc1(i,2)-X_Yr1(j,2))^2<(d_broadcast1(i))^2 &&i~=j
                Vtrue1{i} = [Vtrue1{i} j];
            end
        end
    end
    num_true1 = 0;     %真实干扰边的个数
    for i = 1:Q_true1
        num_true1 = num_true1 + length(Vtrue1{i});
    end
    Q_temp = Q_true;   %不加这句编译就出错!
    Q_temp1 = Q_true1;
    mis_edge = abs(num_true1-num_true);
    pp = mis_edge/num_true;
    data_pp = [data_pp pp]
    slot = slot + M + 1;
    save data_construction,'data_pp';
end

toc;
