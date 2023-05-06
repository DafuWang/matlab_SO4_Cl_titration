function   [V_L,norm_V_L,V_L_dt,norm_V_L_dt,...
            Lambda_dt_1,Lambda_dt_2,delta_Lambda_dt_2_1,delta_theta_2_1,...
                                            norm_Lambda_dt_2_1,norm_delta_theta_2_1]...
            =NO3_Norm_paramrters(k,upsilon,t_1,t_2,V3, ...
                                 Lambda_stage1,Lambda_stage2, ...
                                 Lambda_dt_stage1,Lambda_dt_stage2)     




            Volume_all=1000*upsilon*[t_1,t_2];
            Lambda_all=[Lambda_stage1,Lambda_stage2];
            V_L=[Volume_all',Lambda_all'];

            norm_Volume_all=Volume_all/(V3*1000);
            norm_Lambda_all=(Lambda_all-Lambda_stage1(end))/(abs(Lambda_stage1(1)-Lambda_stage1(end)));
            norm_V_L=[norm_Volume_all',norm_Lambda_all'];

            color=['k','r','g','b','c','m','y'];
            figure(1);
            plot(Volume_all,Lambda_all,'LineWidth',1.5,'Color',color(k),'LineStyle',':')
            set(0,'defaultfigurecolor','w');%设定背景颜色为白色
            set(gca,'FontSize',10,'Fontname', 'Times New Roman');
            xlabel('$Volume\:of\:titration (mL)$','Interpreter','latex','FontSize',15);
            ylabel('$\Lambda (mS/cm)$','Interpreter','latex','FontSize',15);
            set(gcf,'position',[360,198,560,420]);
            set(gca,'position',[0.1,0.1,0.88,0.88]);
            hold on


            N=0;
            Volume_dt_all=1000*upsilon*[t_1(1:end-N),t_2(N+1:end)];
            Lambda_dt_all=[Lambda_dt_stage1(1:end-N),Lambda_dt_stage2(N+1:end)];
            V_L_dt=[Volume_dt_all',Lambda_dt_all'];


            norm_Volume_dt_all=Volume_dt_all/(V3*1000);
            norm_Lambda_dt_all=Lambda_dt_all*V3*1000/(abs(Lambda_stage1(end)-Lambda_stage1(1)));
            norm_V_L_dt=[norm_Volume_dt_all',norm_Lambda_dt_all'];


            
            figure(2)
            plot(Volume_dt_all,Lambda_dt_all,'LineWidth',1.5,'Color',color(k),'LineStyle',':')
            set(0,'defaultfigurecolor','w');%设定背景颜色为白色
            set(gca,'FontSize',10,'Fontname', 'Times New Roman');
            xlabel('$Volume\:of\:titration (mL)$','Interpreter','latex','FontSize',15);
            ylabel('$\dot \Lambda(mS/(cm \cdot mL))$','Interpreter','latex','FontSize',15);
            set(gcf,'position',[360,198,560,420]);
            set(gca,'position',[0.11,0.1,0.87,0.86]);
            hold on


            Lambda_dt_1=Lambda_dt_stage1(end);
            Lambda_dt_2=Lambda_dt_stage2(1);
            delta_Lambda_dt_2_1=Lambda_dt_2-Lambda_dt_1;

            norm_Lambda_dt_1=Lambda_dt_stage1(end)*V3*1000/(abs(Lambda_stage1(end)-Lambda_stage1(1)));
            norm_Lambda_dt_2=Lambda_dt_stage2(1)*V3*1000/(abs(Lambda_stage1(end)-Lambda_stage1(1)));
            norm_Lambda_dt_2_1=norm_Lambda_dt_2-norm_Lambda_dt_1;


            theta_stage1=atan(Lambda_dt_1)*180/pi;%未归一化滴定曲线滴定终点左切线夹角，角度值
            theta_stage2=atan(Lambda_dt_2)*180/pi;%未归一化滴定曲线滴定终点右切线夹角，角度值

            norm_theta_stage1=atan(norm_Lambda_dt_1)*180/pi;%归一化滴定曲线滴定终点左切线夹角，角度值
            norm_theta_stage2=atan(norm_Lambda_dt_2)*180/pi;%归一化滴定曲线滴定终点右切线夹角，角度值


            delta_theta_2_1=theta_stage2-theta_stage1;%未归一化滴定曲线滴定终点左右切线夹角，角度值
            norm_delta_theta_2_1=norm_theta_stage2-norm_theta_stage1;%归一化滴定曲线滴定终点左右切线夹角，角度值         

