function [Lambda,Lambda_dt]=Stage3_AgNO3(bulk,SBNCAN,V,upsilon,t1,t2,t)
% 计算除了硫酸根和滴定离子外其他离子对离子强度的贡献,及其导数,每个时间对应一个数
N=size(bulk.z,2);
n=size(SBNCAN.z,2);
    Gamma=0;Gamma_dt=0;
for i=1:N
    a0=V * bulk.C0(i);
    a0_dt=0;
    b=V + upsilon .* t;
    b_dt=upsilon;

    Gamma=   Gamma    + 0.5 * bulk.z(i)^2. * a0 ./ b;
    Gamma_dt=Gamma_dt + 0.5 * bulk.z(i)^2. *(a0_dt.*b - a0.*b_dt) ./ b.^2;
end

%第1项  代表SO4没有


%第2项 代表Ba
    a2=upsilon .* (t2-t1) .* SBNCAN.C0(2) ;
    a2_dt=0;

    Gamma=   Gamma    + 0.5 *  SBNCAN.z(2)^2 .*a2 ./ b;
    Gamma_dt=Gamma_dt + 0.5 *  SBNCAN.z(2)^2 .*(a2_dt.*b - a2.*b_dt) ./ b.^2;
%第3项 代表来自Ba(NO3)2的NO3
    a3=upsilon .* t2 .* SBNCAN.C0(3);
    a3_dt=0;

    Gamma=   Gamma    + 0.5 *  SBNCAN.z(3)^2 .* a3 ./ b;
    Gamma_dt=Gamma_dt + 0.5 *  SBNCAN.z(3)^2 .* (a3_dt.*b - a3.*b_dt) ./ b.^2;
%第4项 代表Cl
    a4=V * SBNCAN.C0(4)-upsilon .* (t-t2) .* SBNCAN.C0(5);
    a4_dt=-upsilon .* SBNCAN.C0(5);

    Gamma=   Gamma    + 0.5 *SBNCAN.z(4)^2. * a4 ./ b;
    Gamma_dt=Gamma_dt + 0.5 * SBNCAN.z(4)^2. * (a4_dt.*b - a4.*b_dt) ./ b.^2;

% 第5项 代表Ag没有

% 第6项 代表AgNO3的NO3没有
    a6=upsilon .* (t-t2) .* SBNCAN.C0(6);
    a6_dt=upsilon .* SBNCAN.C0(6);

    Gamma=   Gamma    + 0.5 *SBNCAN.z(4)^2. * a6 ./ b;
    Gamma_dt=Gamma_dt + 0.5 * SBNCAN.z(4)^2. * (a6_dt.*b - a6.*b_dt) ./ b.^2;
%% 计算摩尔电导率及其导数，每种离子每个时间对应一个数
for i=1:N
lambda_bulk(i,:)=bulk.lambda_infinity(i) ./ (1+ bulk.G(i) .* sqrt(Gamma));
lambda_bulk_dt(i,:)=-0.5*bulk.lambda_infinity(i).*bulk.G(i).*Gamma_dt./(sqrt(Gamma).*(1+bulk.G(i).*sqrt(Gamma)).^2);
end
%第1项,代表溶液中原有的SO4
i=i+1;j=1;
lambda_SBNCAN(i,1:size(t,2))=0;
lambda_SBNCAN_dt(i,1:size(t,2))=0;

%第2项 代表Ba
i=i+1;j=j+1;
lambda_SBNCAN(i,:)= SBNCAN.lambda_infinity(j) ./ (1+  SBNCAN.G(j) .* sqrt(Gamma));
lambda_SBNCAN_dt(i,:)=-0.5*SBNCAN.lambda_infinity(j).*SBNCAN.G(j).*Gamma_dt./(sqrt(Gamma).*(1+SBNCAN.G(j).*sqrt(Gamma)).^2);
%第3项 代表来自Ba(NO3)2的NO3
i=i+1;j=j+1;
lambda_SBNCAN(i,:)=SBNCAN.lambda_infinity(j) ./ (1+ SBNCAN.G(j) .* sqrt(Gamma));
lambda_SBNCAN_dt(i,:)=-0.5*SBNCAN. lambda_infinity(j).*SBNCAN.G(j).*Gamma_dt./(sqrt(Gamma).*(1+SBNCAN.G(j).*sqrt(Gamma)).^2);

%第4项，代表溶液中原有的Cl
i=i+1;j=j+1;
lambda_SBNCAN(i,:)=SBNCAN.lambda_infinity(j) ./ (1+ SBNCAN.G(j) .* sqrt(Gamma));
lambda_SBNCAN_dt(i,:)=-0.5*SBNCAN.lambda_infinity(j).*SBNCAN.G(j).*Gamma_dt./(sqrt(Gamma).*(1+SBNCAN.G(j).*sqrt(Gamma)).^2);

%第5项，代表来自AgNO3的Ag
i=i+1;j=j+1;
lambda_SBNCAN(i,1:size(t,2))=0;
lambda_SBNCAN_dt(i,1:size(t,2))=0;

%第4项，代表来自AgNO3的NO3
i=i+1;j=j+1;
lambda_SBNCAN(i,:)=SBNCAN.lambda_infinity(j) ./ (1+ SBNCAN.G(j) .* sqrt(Gamma));
lambda_SBNCAN_dt(i,:)=-0.5*SBNCAN.lambda_infinity(j).*SBNCAN.G(j).*Gamma_dt./(sqrt(Gamma).*(1+SBNCAN.G(j).*sqrt(Gamma)).^2);
%% 计算溶液电导率及其导数，每个时间对应一个数
% 计算除了硫酸根和滴定离子外所设计的其他参数，

Lambda=0;Lambda_dt=0;
for i=1:N
    A0=lambda_bulk(i,:) .* V * bulk.C0(i);
    A0_dt=lambda_bulk_dt(i,:) .* V * bulk.C0(i);
    
    Lambda=Lambda + bulk.z(i) .* A0 ./ b;
    Lambda_dt=Lambda_dt+bulk.z(i).*(A0_dt.*b - A0.*b_dt) ./ b.^2;
end
%第1项 没有SO4的贡献
i=i+1;


%第2项 Ba的贡献
i=i+1;
A2=lambda_SBNCAN(i,:) .* a2;
A2_dt=lambda_SBNCAN_dt(i,:).*a2 + lambda_SBNCAN(i,:).*a2_dt;

Lambda=Lambda + SBNCAN.z(1) .*  A2 ./ b;
Lambda_dt=Lambda_dt + SBNCAN.z(1) .* (A2_dt.*b - A2.*b_dt) ./ b.^2;
%第3项参数 NO3的贡献
i=i+1;j=3;
A3=lambda_SBNCAN(i,:) .*a3;
A3_dt= lambda_SBNCAN_dt(i,:).*a3 + lambda_SBNCAN(i,:).*a3_dt;

Lambda=Lambda + SBNCAN.z(j) .* A3 ./ b;
Lambda_dt=Lambda_dt + SBNCAN.z(j) .* (A3_dt.*b - A3.*b_dt) ./ b.^2;

% 第4项 Cl的贡献
i=i+1;j=4;
A4=lambda_SBNCAN(i,:) .* a4;
A4_dt=lambda_SBNCAN_dt(i,:).*a4 + lambda_SBNCAN(i,:).*a4_dt;

Lambda=Lambda + SBNCAN.z(j) .* A4 ./ b;
Lambda_dt=Lambda_dt + SBNCAN.z(j) .* (A4_dt.*b - A4.*b_dt) ./ b.^2;

% 第5项没有
i=i+1;
% 第6项没有
i=i+1;
A6=lambda_SBNCAN(i,:) .* a6;
A6_dt=lambda_SBNCAN_dt(i,:).*a6 + lambda_SBNCAN(i,:).*a6_dt;

Lambda=Lambda + SBNCAN.z(j) .* A6 ./ b;
Lambda_dt=Lambda_dt + SBNCAN.z(j) .* (A6_dt.*b - A6.*b_dt) ./ b.^2;
Lambda_dt=Lambda_dt./upsilon/1000;








