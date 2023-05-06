%%计算电导滴定的滴定液浓度、待测液体积最佳区域图，以及其它离子对最佳滴定区的影响
%Lambda代表溶液电导率，lambda是离子摩尔电导率
%Lambda_dt代表溶液电导率关于t的导数，lambda_dt是离子摩尔电导率关于t的导数

clear
clc
close all

ion=ion_lambda_G_z();%离子有关的电荷和稀溶液摩尔电导率及离子电导率经验系数

%滴定参数设定
upsilon=0.008e-3;%滴定速率
%SO4浓度设定，以Na2SO4的形式
ion.SO4.C0_Na2SO4=0.01;%基准被滴定硫酸钠溶液设定
ion.Na.C0_Na2SO4=2*ion.SO4.C0_Na2SO4;%滴定液中的与硫酸根相应的Na
C0_Na2SO4=ion.SO4.C0_Na2SO4;
%Cl浓度设定，以NaCl的形式
ion.Cl.C0_NaCl=0.01;%基准被滴定硫酸钠溶液设定
ion.Na.C0_NaCl=ion.Cl.C0_NaCl;%滴定液中的与硫酸根相应的Na
C0_NaCl=ion.Cl.C0_NaCl;

% 用于测试SO4的滴定液BaN03_2的设定
C0_BaNO3_2=0.02;%0.005:0.005:0.1;%是一个范围也可以是一个值
% 用于测试Cl的滴定液AgN03的设定
C0_AgNO3=0.05;%0.005:0.005:0.1;

%% 第1个循环是滴定液浓度的变化，第2个循环是其它离子的浓度的变化，第3个循环是待测液体积的变化
i=0;
for C_BaNO3_2=C0_BaNO3_2%
    i=i+1;
    ion.Ba.C0_BaNO3_2=C_BaNO3_2;%滴定浓度
    ion.NO3.C0_BaNO3_2=2*ion.Ba.C0_BaNO3_2;%滴定浓度
    
    ion.Ag.C0_AgNO3=C0_AgNO3;%Ag滴定浓度
    ion.NO3.C0_AgNO3=C0_AgNO3;%NO3滴定浓度
    %添加进去的影响离子浓度设定
    j=0;
    C0_bulk_NaNO3=0;%0.5:0.5:0.5;
    for ion_bulk_X=C0_bulk_NaNO3
        j=j+1;
        ion.K.C0_bulk=0;
        ion.Ca.C0_bulk=0;
        ion.NO3.C0_bulk=ion_bulk_X;%溶液中原有的NO3

        ion.H.C0_bulk=0;
        ion.Na.C0_bulk1=ion_bulk_X;%溶液中不是来自Na2SO4和NaCl
        ion.OH.C0_bulk=0;

        ion.Na.C0_bulk=ion.Na.C0_bulk1+2*C0_Na2SO4+C0_NaCl;%待测液中总的Na
  
        [bulk,SBNCAN]=ion_z_C0_G_lambda_infinity_bulk(ion);
        k=0;
        V_variation=5e-3:5e-3:30e-3;
        for V=V_variation%改变滴定溶液体积
            k=k+1;
%% SO4滴定      
            V1_SO4=V*C0_Na2SO4/C_BaNO3_2;%理论滴定体积
            %% SO4滴定第1阶段
            t1_SO4=round(V1_SO4/upsilon);%第1阶段结束的时间
            t_1_SO4=0:t1_SO4;%第1阶段的时间
            [Lambda_stage1_SO4,Lambda_dt_stage1_SO4,lambda_SO4_t0]=Stage1_BaNO3_2(bulk,SBNCAN,V,upsilon,t_1_SO4);%lambda_SO4_t0产生的电导率

            %% SO4滴定第2阶段
             t2_SO4=round(1.1*V1_SO4/upsilon);%假设第2阶段结束的时间
             t_2_SO4=t1_SO4:t2_SO4;%第2阶段的时间

            [Lambda_stage2_SO4,Lambda_dt_stage2_SO4]=Stage2_BaNO3_2(bulk,SBNCAN,V,upsilon,t1_SO4, t_2_SO4);
          
            %% 数据整合
            [V_L_SO4{i,k},norm_V_L_SO4{i,j,k},V_L_dt_SO4{i,k},norm_V_L_dt_SO4{i,j,k},...
             Lambda_dt_1_SO4(i,j,k),Lambda_dt_2_SO4(i,j,k),delta_Lambda_dt_2_1_SO4(i,j,k),delta_theta_2_1_SO4(i,j,k),...
                                                           norm_Lambda_dt_2_1_SO4(i,j,k),norm_delta_theta_2_1_SO4(i,j,k)]...
            =SO4_Norm_paramrters(k,upsilon,t_1_SO4,t_2_SO4,V1_SO4,...
                                 Lambda_stage1_SO4,Lambda_stage2_SO4, ...
                                 Lambda_dt_stage1_SO4,Lambda_dt_stage2_SO4);

%% NO3滴定      
            V1_NO3=V*C0_NaCl/C0_AgNO3;%理论滴定体积
            %% NO3滴定第1阶段
            t1_NO3=t2_SO4+round(V1_NO3/upsilon);%第3阶段结束的时间
            t_1_NO3=t2_SO4:t1_NO3;%第3阶段的时间

            [Lambda_stage1_NO3,Lambda_dt_stage1_NO3]=Stage3_AgNO3(bulk,SBNCAN,V,upsilon,t1_SO4,t2_SO4,t_1_NO3);%lambda_SO4_t0产生的电导率
            %% NO3滴定第2阶段
            t2_NO3=t2_SO4+round(1.5*V1_NO3/upsilon);%第3阶段结束的时间
            t_2_NO3=t1_NO3:t2_NO3;%第3阶段的时间
            [Lambda_stage2_NO3,Lambda_dt_stage2_NO3]=Stage4_AgNO3(bulk,SBNCAN,V,upsilon,t1_SO4,t2_SO4,t1_NO3,t_2_NO3);%lambda_SO4_t0产生的电导率
            %% 数据整合
            [V_L_NO3{i,k},norm_V_L_NO3{i,j,k},V_L_dt_NO3{i,k},norm_V_L_dt_NO3{i,j,k},...
             Lambda_dt_1_NO3(i,j,k),Lambda_dt_2_NO3(i,j,k),delta_Lambda_dt_2_1_NO3(i,j,k),delta_theta_2_1_NO3(i,j,k),...
                                                           norm_Lambda_dt_2_1_NO3(i,j,k),norm_delta_theta_2_1_NO3(i,j,k)]...
            =NO3_Norm_paramrters(k,upsilon,t_1_NO3,t_2_NO3,V1_NO3,...
                                 Lambda_stage1_NO3,Lambda_stage2_NO3, ...
                                 Lambda_dt_stage1_NO3,Lambda_dt_stage2_NO3);




        end
    end
end


 
% for i=1:size(C0_Ba2N03,2)
%     V_x=V_variation'*1000;%待测液体积的变化
% 
%     figure(7)%非归一化滴定曲线滴定终点左右切线的夹角随体积的变化
%     y1(:,i)=delta_theta_5_3_SO4(i,j,:);
%     plot(V_x,y1(:,i))
%     hold on
% 
%     figure(8)%归一化滴定曲线滴定终点左右切线的夹角随体积的变化
%     y2(:,i)=norm_delta_theta_5_3_SO4(i,j,:);
%     plot(V_x,y2(:,i))
%     hold on
% 
% end

table0=[Lambda_dt_1_SO4,Lambda_dt_2_SO4,delta_Lambda_dt_2_1_SO4,delta_theta_2_1_SO4,norm_Lambda_dt_2_1_SO4,norm_delta_theta_2_1_SO4];
table1(:,:)=table0(1,:,:);
table_SO4=table1';


table0_NO3=[Lambda_dt_1_NO3,Lambda_dt_2_NO3,delta_Lambda_dt_2_1_NO3,delta_theta_2_1_NO3,norm_Lambda_dt_2_1_NO3,norm_delta_theta_2_1_NO3];
table1_NO3(:,:)=table0_NO3(1,:,:);
table_NO3=table1_NO3';














