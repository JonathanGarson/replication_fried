%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: read_c.m
% Bay: Stephie Fried
% Date: Summer 2021
% Purpose: Reads and processes C results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all

datapath = strcat('C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\c_programs\results_1109');
bpath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output';
cd(bpath);

%%
% Baseline results
cd(datapath);

%Flags
%b- baseline
%np - no fema aid and no subsidy
%ns - fema aid, no subsidy 
%nr - no risk 
%ccna - climate change, no adaptation
%cc - climate change 
%ccnr- climate change, no risk
%hr-  higher rho, used to compute the elasticities 

exper_vec1 = {'b', 'np', 'ns', 'nr', 'hr'};
exper_vec2 = {'ccna', 'cc', 'ccnr'};
nc = 5; %number of climate change experiments minus 1

%Import parameter values
M = importdata('params.dat');
deltak= M(1); zeta= M(2); beta= M(3); theta= M(4); Omega_hh= M(5); Omegak_Omegah = M(6);
sigma= M(7); psiky= M(8); psiho = M(9); rho_l= M(10); rho_h = M(11);  eta= M(12); lambda= M(13);
deltah= M(14); r= M(15); tol = M(16); tfp_l = M(17); tfp_h = M(18);
alpha = M(19);  dtol = M(20); gtol  = M(21); ftol = M(22);  phi_r = M(23); hpmin = M(24);
psihr = M(25); Ah = M(26); deltaa = M(27); chi = M(28); psikr = M(29);

%Import moments and targets
M = importdata('targets_b.dat');
wealth_output = M(1); aid_output= M(2); subs_output = M(3); housing_capital = M(4); owned_total = M(5);
rel_gdp = M(6); insurance_house_damage = M(7); housing_damage_total_damage = M(8);
owner_frac = M(9); fema_damage= M(10); renter_aid_total_aid= M(11); owner_aid_total_aid= M(12);
rental_housing_aid_total_aid= M(13); rel_fema_aid_owners= M(14); rel_fema_aid_renters= M(15);
rel_fema_aid_ffirms= M(16); wealth_output_target= M(17); aid_output_target= M(18);
subs_output_target= M(19); housing_capital_target= M(20); owned_total_target= M(21); rel_gdp_target= M(22);
insurance_house_damage_target= M(23); housing_damage_total_damage_target= M(24); owner_frac_target= M(25);
rel_damage_target= M(26); renter_aid_total_aid_target= M(27); owner_aid_total_aid_target= M(28);
rental_housing_aid_total_aid_target= M(29); fema_damage_target= M(30); rel_damage_renters_target= M(31);
rel_damage_firms_target= M(32); owner_renter_inc_target= M(33); owner_renter_inc= M(34);
i1= M(35); i2= M(36); i3= M(37); i4= M(38); i5 = M(39); i1_target= M(40); i2_target= M(41); i3_target= M(42);
i4_target= M(43); i5_target = M(44);
w1 = M(45);w2= M(46);w3= M(47); w4 = M(48);w5 = M(49);w1_target= M(50); w2_target = M(51);w3_target= M(52);
w4_target = M(53);w5_target= M(54); rel_aid_target = M(55); rel_aid = M(56);

% Create grids
M = importdata('grids.dat');
nhs = M(1);  nm = M(2); nz= M(3); nq =2;
hsmax= M(4); hsmin= M(5); mmax= M(6); mmin= M(7);

mgrid = importdata('mgrid.dat');
hsgrid = importdata('hsgrid.dat');
zgrid = importdata('zgrid.dat');
qgrid = [1,2];
tmat_vec = importdata('tmat.dat');

counter =0;
for i = 1:nz
    for  j = 1:nz
        counter =counter +1;
        tmat(i,j) = tmat_vec(counter);
    end
end

%Find exp1ected z next period, given qi and zi
Ez = zeros(nz, 1);
temp = zeros(nz,1);
for zi = 1:nz
    for z2i =1:nz
        Ez(zi) = Ez(zi) + tmat(zi, z2i)*zgrid(z2i);
        temp (zi) = temp(zi) + tmat(zi, z2i); %check that temp equals 1
    end
end

n = nhs*nm*nz*nq;

for i = 1:n
    j = i-1;
    qind(i) = floor(j/(nm*nhs*nz)) +1;
    zind(i)  = mod(floor(j/(nm*nhs)), nz)+1;
    hsind(i)  = mod(floor(j/(nm)), nhs)+1;
    mind(i) = mod(j, nm)+1;
end


%Save baseline values of parameters and targets
cd(bpath);
save params_targs_base;

reg_vec1 = {'l', 'h'};

%Import baseline and non-climate change exp1eriment results
for k = 1:length(exper_vec1)
    k
    exper = exper_vec1{k};
    clearvars -regexp _vec$
    
    %Import value functions
    cd(datapath);
    v_l_vec = importdata(strcat('v0_', exper, '0', '.dat'));
    f_l_vec = importdata(strcat('f1_', exper, '0', '.dat'));
    v_h_vec = importdata(strcat('v0_', exper, '1', '.dat'));
    f_h_vec = importdata(strcat('f1_', exper, '1', '.dat'));
    
    %Find welfare
    wel = 0.5*f_l_vec'*v_l_vec + 0.5*f_h_vec'*v_h_vec;
    wel_l =f_l_vec'*v_l_vec;
    wel_h =f_h_vec'*v_h_vec;
    eval(strcat('wel_', exper, '= wel;'));
    eval(strcat('wel_l_', exper, '= wel_l;'));
    eval(strcat('wel_h_', exper, '= wel_h;'));
    
    for m = 1:2
        reg = reg_vec1{m};
        eval(strcat('v_', reg, '_vec_', exper, '= v_', reg, '_vec;'));
        eval(strcat('f_', reg, '_vec_', exper, '= f_', reg, '_vec;'));
    end
    
    %Import policy functions
    hppol_l_vec = importdata(strcat('hppol_', exper, '0', '.dat'));
    apol_l_vec = importdata(strcat('apol_', exper, '0', '.dat'));
    hsrpol_l_vec = importdata(strcat('hsrpol_', exper, '0', '.dat'));
    hppol_h_vec = importdata(strcat('hppol_', exper, '1', '.dat'));
    apol_h_vec = importdata(strcat('apol_', exper, '1', '.dat'));
    hsrpol_h_vec = importdata(strcat('hsrpol_', exper, '1', '.dat'));
    hapol_l_vec = importdata(strcat('hapol_', exper, '0', '.dat'));
    xpol_l_vec = importdata(strcat('xpol_', exper, '0', '.dat'));
    hapol_h_vec = importdata(strcat('hapol_', exper, '1', '.dat'));
    xpol_h_vec = importdata(strcat('xpol_', exper, '1', '.dat'));
    cd(bpath);
    
    for m = 1:2
        reg = reg_vec1{m};
        eval(strcat('hppol_', reg, '_vec_', exper, '= hppol_', reg, '_vec;'));
        eval(strcat('hapol_', reg, '_vec_', exper, '= hapol_', reg, '_vec;'));
        eval(strcat('apol_', reg, '_vec_', exper, '= apol_', reg, '_vec;'));
        eval(strcat('xpol_', reg, '_vec_', exper, '= xpol_', reg, '_vec;'));
        eval(strcat('hsrpol_', reg, '_vec_', exper, '= hsrpol_', reg, '_vec;'));
    end
    
    %Unwind value and policy functions
    for m = 1:2
        reg = reg_vec1{m};
        eval(strcat('v_vec = v_', reg, '_vec;'));
        eval(strcat('f_vec = f_', reg, '_vec;'));
        eval(strcat('hppol_vec = hppol_', reg, '_vec;'));
        eval(strcat('hapol_vec = hapol_', reg, '_vec;'));
        eval(strcat('apol_vec = apol_', reg, '_vec;'));
        eval(strcat('xpol_vec = xpol_', reg, '_vec;'));
        eval(strcat('hsrpol_vec = hsrpol_', reg, '_vec;'));
        
        for i = 1:n
            v(qind(i), zind(i), hsind(i), mind(i)) = v_vec(i);
            f(qind(i), zind(i), hsind(i), mind(i)) = f_vec(i);
            hppol(qind(i), zind(i), hsind(i), mind(i)) = hppol_vec(i);
            hapol(qind(i), zind(i), hsind(i), mind(i)) = hapol_vec(i);
            apol(qind(i), zind(i), hsind(i), mind(i)) = apol_vec(i);
            xpol(qind(i), zind(i), hsind(i), mind(i)) = xpol_vec(i);
            hsrpol(qind(i), zind(i), hsind(i), mind(i)) = hsrpol_vec(i);
            
        end
        
        eval(strcat('v_', reg, '_', exper, '= v;'));
        eval(strcat('f_', reg, '_', exper, '= f;'));
        eval(strcat('hppol_', reg, '_', exper, '= hppol;'));
        eval(strcat('apol_', reg, '_', exper, '= apol;'));
        eval(strcat('xpol_', reg, '_', exper, '= xpol;'));
        eval(strcat('hsrpol_', reg, '_', exper, '= hsrpol;'));
        eval(strcat('hapol_', reg, '_', exper, '= hapol;'));
    end
    
end

%Import climate change exp1eriment results
for k = 1:length(exper_vec2)
    k
    exper = exper_vec2{k};
    for si =0:1:nc+1 %Read in each climate change scenario
        clearvars -regexp _vec$
        if (si == nc+1) si = 9; end
        if ( (exper =="cc" && si ==9)  || (exper =="ccnr" && si ==9))
            5+5;
        else
            
            %Import value functions
            cd(datapath);
            v_l_vec = importdata(strcat('v0_', exper, '0', '_', num2str(si), '.dat'));
            f_l_vec = importdata(strcat('f1_', exper, '0',  '_', num2str(si), '.dat'));
            v_h_vec = importdata(strcat('v0_', exper, '1',  '_', num2str(si), '.dat'));
            f_h_vec = importdata(strcat('f1_', exper, '1',  '_', num2str(si), '.dat'));
            
            %Find welfare
            wel = 0.5*f_l_vec'*v_l_vec + 0.5*f_h_vec'*v_h_vec;
            wel_l =f_l_vec'*v_l_vec;
            wel_h =f_h_vec'*v_h_vec;
            eval(strcat('wel_', exper,num2str(si), '= wel;'));
            eval(strcat('wel_l_', exper,num2str(si), '= wel_l;'));
            eval(strcat('wel_h_', exper,num2str(si), '= wel_h;'));
            
            for m = 1:2
                reg = reg_vec1{m};
                eval(strcat('v_', reg, '_vec_', exper, num2str(si), '= v_', reg, '_vec;'));
                eval(strcat('f_', reg, '_vec_', exper, num2str(si), '= f_', reg, '_vec;'));
            end
            
            
            %Import policy functions
            hppol_l_vec = importdata(strcat('hppol_', exper, '0', '_', num2str(si), '.dat'));
            apol_l_vec = importdata(strcat('apol_', exper, '0', '_', num2str(si), '.dat'));
            hsrpol_l_vec = importdata(strcat('hsrpol_', exper, '0', '_', num2str(si), '.dat'));
            hppol_h_vec = importdata(strcat('hppol_', exper, '1', '_', num2str(si), '.dat'));
            apol_h_vec = importdata(strcat('apol_', exper, '1', '_', num2str(si), '.dat'));
            hsrpol_h_vec = importdata(strcat('hsrpol_', exper, '1', '_', num2str(si), '.dat'));
            hapol_l_vec = importdata(strcat('hapol_', exper, '0', '_', num2str(si), '.dat'));
            xpol_l_vec = importdata(strcat('xpol_', exper, '0', '_', num2str(si), '.dat'));
            hapol_h_vec = importdata(strcat('hapol_', exper, '1', '_',num2str(si), '.dat'));
            xpol_h_vec = importdata(strcat('xpol_', exper, '1', '_', num2str(si), '.dat'));
            cd(bpath);
            
            for m = 1:2
                reg = reg_vec1{m};
                eval(strcat('hppol_', reg, '_vec_', exper, num2str(si),  '= hppol_', reg, '_vec;'));
                eval(strcat('hapol_', reg, '_vec_', exper, num2str(si), '= hapol_', reg, '_vec;'));
                eval(strcat('apol_', reg, '_vec_', exper, num2str(si), '= apol_', reg, '_vec;'));
                eval(strcat('xpol_', reg, '_vec_', exper, num2str(si), '= xpol_', reg, '_vec;'));
                eval(strcat('hsrpol_', reg, '_vec_', exper, num2str(si), '= hsrpol_', reg, '_vec;'));
                
            end
            
            %Unwind value and policy functions
            for m = 1:2
                reg = reg_vec1{m};
                eval(strcat('v_vec = v_', reg, '_vec;'));
                eval(strcat('f_vec = f_', reg, '_vec;'));
                eval(strcat('hppol_vec = hppol_', reg, '_vec;'));
                eval(strcat('hapol_vec = hapol_', reg, '_vec;'));
                eval(strcat('apol_vec = apol_', reg, '_vec;'));
                eval(strcat('xpol_vec = xpol_', reg, '_vec;'));
                eval(strcat('hsrpol_vec = hsrpol_', reg, '_vec;'));
                
                
                for i = 1:n
                    v(qind(i), zind(i), hsind(i), mind(i)) = v_vec(i);
                    f(qind(i), zind(i), hsind(i), mind(i)) = f_vec(i);
                    hppol(qind(i), zind(i), hsind(i), mind(i)) = hppol_vec(i);
                    hapol(qind(i), zind(i), hsind(i), mind(i)) = hapol_vec(i);
                    apol(qind(i), zind(i), hsind(i), mind(i)) = apol_vec(i);
                    xpol(qind(i), zind(i), hsind(i), mind(i)) = xpol_vec(i);
                    hsrpol(qind(i), zind(i), hsind(i), mind(i)) = hsrpol_vec(i);
                end
                
                eval(strcat('v_', reg, '_', exper,num2str(si),  '= v;'));
                eval(strcat('f_', reg, '_', exper, num2str(si), '= f;'));
                eval(strcat('hppol_', reg, '_', exper,num2str(si),  '= hppol;'));
                eval(strcat('apol_', reg, '_', exper,num2str(si),  '= apol;'));
                eval(strcat('xpol_', reg, '_', exper,num2str(si),  '= xpol;'));
                eval(strcat('hsrpol_', reg, '_', exper,num2str(si),  '= hsrpol;'));
                eval(strcat('hapol_', reg, '_', exper,num2str(si),  '= hapol;'));
                
            end
        end
    end
end

save value_policy_funcs

%% Aggregates

vars = {'Kp_reg', 'Ka_reg', 'Kptilde_reg', 'Dk_reg', 'w_reg','Hp_reg', 'Hrp_reg', 'Ha_reg', 'A_reg', ...
    'Hra_reg', 'Ytilde_reg', 'Kp', 'Hp', 'Hrp', 'frac_adapt', 'mdist_sim_peeps', 'Fomega_own', 'renters_reg',...
    'frac_adapt_k', 'frac_adapt_h', 'Ytilde', 'D', 'frac_adapt_reg', 'frac_adapt_k_reg', ...
    'frac_adapt_h_reg', 'tau', 'Ha', 'A', 'tdiff', 'Ka', 'Omegah_reg', 'theta', 'psiky', 'psikr', ...
    'eta','rho_reg', 'Dhr_reg', 'ph_reg', 'Omegak_reg', 'psiho', 'GDP', 'GDP_reg', 'Dhr', 'Dhp', 'Dk',...
    'Dhp_reg', 'total_cap', 'ar_reg', 'ay_reg', 'Y_reg', 'D_reg', 'D', 'D_output', 'D_cap', 'x', 'total_cap_reg',...
    'D_totalK_reg',  'aid_reg', 'subs_reg', 'Hra', 'Hrp', 'psihr', 'owners_reg', 'iter', 'mdist', 'iter_sim_peeps', ...
    'f_vec', 'hp_vec', 'ha_vec', 'a_vec', 'v_vec', 'x_vec', 'Omegah_reg', 'f_ho_vec', 'ah_vec', 'dhp_vec',...
    'ah', 'ah_reg', 'x_dhp_vec', 'wealth', 'wealth_vec', 'x_dhp', 'x_dhp_reg',  'aid', 'subs', 'renters', 'wealth_reg'};


vars2 = {'hppol_l_vec', 'hppol_h_vec', 'hapol_l_vec', 'hapol_h_vec', 'apol_l_vec', 'apol_h_vec', 'v_l_vec', 'v_h_vec',...
    'f_l_vec', 'f_h_vec', 'xpol_l_vec', 'xpol_h_vec'};

%Non Climate change exp1eriments
for k = 1:length(exper_vec1)
    k
    exper = exper_vec1{k};
    
    for m = 1:2
        reg = reg_vec1{m};
        cd(datapath);
        M = importdata(strcat('agg_', exper, num2str(m-1), '.dat'));
        cd(bpath);
        Kp_reg(m) = M(1); Ka_reg(m) = M(2); Kptilde_reg(m) = M(3); Dk_reg(m) = M(4);  w_reg(m) = M(5);
        Hp_reg(m) = M(6); Hrp_reg(m) = M(7); Ha_reg(m) = M(8); A_reg(m) = M(9);
        x_dhp_reg(m) = M(10); owners_reg(m) = M(11); Hra_reg(m) = M(12); tau = M(13); D_reg(m) = M(14);
        Ika_reg(m) = M(15); theta =M(16); Omegah_reg(m) = M(17); psiky = M(18);
        eta = M(19); tdiff = M(20); ph_reg(m) = M(21); rho_reg(m) = M(22); Dhr_reg(m) = M(23); psikr = M(24);
        ar_reg(m) = M(25); ay_reg(m) = M(26); Omegak_reg(m) = M(27); Dhp_reg(m) = M(28);
        Ytilde_reg(m) = M(29); GDP_reg(m) = M(30); psiho = M(31); psihr = M(32); iter = M(33); mdist = M(34);
        iter_sim_peeps = M(35); mdist_sim_peeps = M(36);
    end
    
    Kp =mean(Kp_reg); Ka = mean(Ka_reg);
    Hp = mean(Hp_reg); Ha = mean(Ha_reg);
    Hrp = mean(Hrp_reg); Hra = mean(Hra_reg);
    A = mean(A_reg);
    Ika = mean(Ika_reg);
    Ihra = deltaa*Hra; Iha = deltaa*Ha;
    subs_reg= eta*deltaa*(Ka_reg + Hra_reg + Ha_reg);
    subs = mean(subs_reg);
    wealth_reg = Hp_reg + Ha_reg./(1+eta) + A_reg;
    
    frac_adapt = (Ka + Hra + Ha)/(Ka +Kp +Hra + Hrp +Ha + Hp);
    frac_adapt_k = Ka/(Ka + Kp);
    frac_adapt_h = (Hra + Ha)/(Hra + Hrp + Ha + Hp);
    GDP = mean(GDP_reg);
    Ytilde = mean(Ytilde_reg);
    
    D_cap = mean(Dk_reg) + mean(Dhp_reg) + mean(Dhr_reg);
    Y_reg = Kp_reg.^alpha;
    D_output_reg = Y_reg - Ytilde_reg;
    D_output = mean(D_output_reg);  %Damage to output
    D_reg = D_output_reg + Dk_reg + Dhp_reg + Dhr_reg + ph_reg.*Ah.*(Dhr_reg + Dhp_reg);
    D = mean(D_reg);
    Dhp = mean(Dhp_reg); Dhr = mean(Dhr_reg); Dk = mean(Dk_reg);
    aid_reg = (psiky*(Dhr_reg + Dk_reg) + psihr*Dhr_reg + psiho*Dhp_reg);
    aid = mean(aid_reg);
    total_cap = Ka + Kp + Hra + Hrp + Ha + Hp;
    total_cap_reg = Ka_reg+ Kp_reg + Hra_reg + Ha_reg + Hp_reg;
    D_totalK_reg = Dk_reg + Dhp_reg + Dhr_reg;
    
    %Regional aggregates
    for i = 1:2
        frac_adapt_reg(i) =  (Ka_reg(i) + Hra_reg(i) + Ha_reg(i))...
            /(Ka_reg(i) +Kp_reg(i) +Hra_reg(i) + Hrp_reg(i) + Ha_reg(i) + Hp_reg(i));
        frac_adapt_k_reg(i) = Ka_reg(i)/(Ka_reg(i) + Kp_reg(i));
        frac_adapt_h_reg(i) = (Hra_reg(i) + Ha_reg(i))...
            /(Hra_reg(i) + Hrp_reg(i) + Ha_reg(i) + Hp_reg(i));
    end
    
    %Homeowner aggregates
    for i = 1:length(vars2)
        eval(strcat(vars2{i}, '=', vars2{i}, '_', exper, ';'));
    end
    
    f_vec = [0.5*f_l_vec; 0.5*f_h_vec];
    hp_vec = [hppol_l_vec; hppol_h_vec];
    ha_vec = [hapol_l_vec; hapol_h_vec];
    a_vec = [apol_l_vec; apol_h_vec];
    v_vec = [v_l_vec; v_h_vec];
    x_vec = [xpol_l_vec; xpol_h_vec];
    Omegah_vec =[Omegah_reg(1)*ones(n,1); Omegah_reg(2)*ones(n,1)];
    
    %Values just for homeowners
    f_ho_vec = f_vec;
    f_ho_vec(hp_vec ==0)=0;
    ah_vec = ha_vec./hp_vec;
    ah_vec(isnan(ah_vec))=0;
    dhp_vec = Omegah_vec.*exp( -theta*(ah_vec).^(1/theta)).*hp_vec;
    x_dhp_vec = x_vec./dhp_vec;
    x_dhp_vec(isnan(x_dhp_vec))=0;
    wealth_vec = hp_vec +ha_vec/(1+eta) + a_vec;
    
    ah = f_ho_vec'*ah_vec./(sum(f_ho_vec));
    wealth = f_ho_vec'*wealth_vec./(sum(f_ho_vec));
    ah_reg(1) = f_ho_vec(1:n)'*ah_vec(1:n)/sum(f_ho_vec(1:n));
    ah_reg(2) = f_ho_vec(n+1:2*n)'*ah_vec(n+1:2*n)/sum(f_ho_vec(n+1:2*n));
    x_dhp = f_ho_vec'*x_dhp_vec./(sum(f_ho_vec));
    x_dhp_reg(1) = f_ho_vec(1:n)'*x_dhp_vec(1:n)/sum(f_ho_vec(1:n));
    x_dhp_reg(2) = f_ho_vec(n+1:2*n)'*x_dhp_vec(n+1:2*n)/sum(f_ho_vec(n+1:2*n));
    x = f_ho_vec'*x_vec/sum(f_ho_vec(1:n));
    Fomega_own = f_ho_vec'*(exp(-theta*ah_vec.^(1/theta)).*Omegah_vec)/sum(f_ho_vec);
    
    renters = 1 - sum(f_ho_vec);
    renters_reg(1) = 1 -2*sum(f_ho_vec(1:n/2));
    renters_reg(2) = 1 -2*sum(f_ho_vec(n/2+1:end));
    
    
    filename = strcat('agg_', exper);
    for i = 1:length(vars)
        str = strcat(vars{i}, '_', exper, '=', vars{i}, ';');
        eval(str);
        str2 =  strcat(vars{i}, '_', exper);
        if i ==1 save(filename, str2);
        else
            save(filename, str2, '-append');
        end
    end
end

%Climate change exp1eriments
for k = 1:length(exper_vec2)
    k
    exper = exper_vec2{k};
    
    for si = 0:1:nc+1
        if (si == nc+1) si = 9; end
        if ( (exper =="cc" && si ==9)  || (exper =="ccnr" && si ==9))
            5+5;
        else
            
            for m = 1:2
                reg = reg_vec1{m};
                cd(datapath);
                M = importdata(strcat('agg_', exper,  num2str(m-1), '_', num2str(si), '.dat'));
                cd(bpath);
                Kp_reg(m) = M(1); Ka_reg(m) = M(2); Kptilde_reg(m) = M(3); Dk_reg(m) = M(4);  w_reg(m) = M(5);
                Hp_reg(m) = M(6); Hrp_reg(m) = M(7); Ha_reg(m) = M(8); A_reg(m) = M(9);
                x_dhp_reg(m) = M(10); owners_reg(m) = M(11); Hra_reg(m) = M(12); tau = M(13); D_reg(m) = M(14);
                Ika_reg(m) = M(15); theta =M(16); Omegah_reg(m) = M(17); psiky = M(18);
                eta = M(19); tdiff = M(20); ph_reg(m) = M(21); rho_reg(m) = M(22); Dhr_reg(m) = M(23); psikr = M(24);
                ar_reg(m) = M(25); ay_reg(m) = M(26); Omegak_reg(m) = M(27); Dhp_reg(m) = M(28);
                Ytilde_reg(m) = M(29); GDP_reg(m) = M(30); psiho = M(31); psihr = M(32); iter = M(33); mdist = M(34);
                iter_sim_peeps = M(35); mdist_sim_peeps = M(36);
            end
            
            Kp =mean(Kp_reg); Ka = mean(Ka_reg);
            Hp = mean(Hp_reg); Ha = mean(Ha_reg);
            Hrp = mean(Hrp_reg); Hra = mean(Hra_reg);
            A = mean(A_reg);
            Ika = mean(Ika_reg);
            Ihra = deltaa*Hra;
            Iha = deltaa*Ha;
            subs_reg= eta*deltaa*(Ka_reg + Hra_reg + Ha_reg);
            wealth_reg = Hp_reg + Ha_reg./(1+eta) + A_reg;
            
            frac_adapt = (Ka + Hra + Ha)/(Ka +Kp +Hra + Hrp +Ha + Hp);
            frac_adapt_k = Ka/(Ka + Kp);
            frac_adapt_h = (Hra + Ha)/(Hra + Hrp + Ha + Hp);
            GDP = mean(GDP_reg);
            Ytilde = mean(Ytilde_reg);
            X_dhp  = mean(x_dhp_reg);
            
            
            %calculate damage
            D_cap = mean(Dk_reg) + mean(Dhp_reg) + mean(Dhr_reg);
            Y_reg = Kp_reg.^alpha;
            D_output_reg = Y_reg - Ytilde_reg;
            D_output = mean(D_output_reg);
            D_reg = D_output_reg + Dk_reg + Dhp_reg + Dhr_reg + ph_reg.*Ah.*(Dhr_reg + Dhp_reg);
            D = mean(D_reg);
            Dhp = mean(Dhp_reg); Dhr = mean(Dhr_reg); Dk = mean(Dk_reg);
            aid_reg = (psiky*(Dhr_reg + Dk_reg) + psihr*Dhr_reg + psiho*Dhp_reg);
            aid = mean(aid_reg);
            total_cap = Ka + Kp + Hra + Hrp + Ha + Hp;
            total_cap_reg = Ka_reg+ Kp_reg + Hra_reg + Ha_reg + Hp_reg;
            D_totalK_reg = Dk_reg + Dhp_reg + Dhr_reg;
            
            %Regional aggregates
            for i = 1:2
                frac_adapt_reg(i) =  (Ka_reg(i) + Hra_reg(i) + Ha_reg(i))...
                    /(Ka_reg(i) +Kp_reg(i) +Hra_reg(i) + Hrp_reg(i) + Ha_reg(i) + Hp_reg(i));
                frac_adapt_k_reg(i) = Ka_reg(i)/(Ka_reg(i) + Kp_reg(i));
                frac_adapt_h_reg(i) = (Hra_reg(i) + Ha_reg(i))...
                    /(Hra_reg(i) + Hrp_reg(i) + Ha_reg(i) + Hp_reg(i));
            end
            
            %Homeowner aggregates
            for i = 1:length(vars2)
                eval(strcat(vars2{i}, '=', vars2{i}, '_', exper, num2str(si), ';'));
            end
            
            f_vec = [0.5*f_l_vec; 0.5*f_h_vec];
            hp_vec = [hppol_l_vec; hppol_h_vec];
            ha_vec = [hapol_l_vec; hapol_h_vec];
            a_vec = [apol_l_vec; apol_h_vec];
            v_vec = [v_l_vec; v_h_vec];
            x_vec = [xpol_l_vec; xpol_h_vec];
            wealth_vec = hp_vec +ha_vec/(1+eta) + a_vec;
            Omegah_vec =[Omegah_reg(1)*ones(n,1); Omegah_reg(2)*ones(n,1)];
            
            %Values just for homeowners
            f_ho_vec = f_vec;
            f_ho_vec(hp_vec ==0)=0;
            ah_vec = ha_vec./hp_vec;
            ah_vec(isnan(ah_vec))=0;
            dhp_vec = Omegah_vec.*exp( -theta*(ah_vec).^(1/theta)).*hp_vec;
            x_dhp_vec = x_vec./dhp_vec;
            x_dhp_vec(isnan(x_dhp_vec))=0;
            
            ah = f_ho_vec'*ah_vec./(sum(f_ho_vec));
            wealth = f_ho_vec'*wealth_vec./(sum(f_ho_vec));
            ah_reg(1) = f_ho_vec(1:n)'*ah_vec(1:n)/sum(f_ho_vec(1:n));
            ah_reg(2) = f_ho_vec(n+1:2*n)'*ah_vec(n+1:2*n)/sum(f_ho_vec(n+1:2*n));
            x_dhp = f_ho_vec'*x_dhp_vec./(sum(f_ho_vec));
            x_dhp_reg(1) = f_ho_vec(1:n)'*x_dhp_vec(1:n)/sum(f_ho_vec(1:n));
            x_dhp_reg(2) = f_ho_vec(n+1:2*n)'*x_dhp_vec(n+1:2*n)/sum(f_ho_vec(n+1:2*n));
            x_dhp_reg(1) = f_ho_vec(1:n)'*x_dhp_vec(1:n)/sum(f_ho_vec(1:n));
            x = f_ho_vec'*x_vec/sum(f_ho_vec(1:n));
            Fomega_own = f_ho_vec'*(exp(-theta*ah_vec.^(1/theta)).*Omegah_vec)/sum(f_ho_vec);
            
            renters = 1 - sum(f_ho_vec);
            renters_reg(1) = 1 -2*sum(f_ho_vec(1:n/2));
            renters_reg(2) = 1 -2*sum(f_ho_vec(n/2+1:end));
            
            
            filename = strcat('agg_', exper, num2str(si));
            for i = 1:length(vars)
                str = strcat(vars{i}, '_', exper, num2str(si), '=', vars{i}, ';');
                eval(str);
                str2 =  strcat(vars{i}, '_', exper, num2str(si));
                if i ==1 save(filename, str2);
                else
                    save(filename, str2, '-append');
                end
            end
        end
    end
end

%% Aggregates by quintile: baseline
load agg_b
load params_targs_base
x2 = 0.05; %threshold for barely unconstrained homeowners
reg_vec =[ones(n, 1); 2*ones(n,1)];
m_vec = [mgrid(mind); mgrid(mind)];
hs_vec = [hsgrid(hsind); hsgrid(hsind)];
z_vec= [zgrid(zind); zgrid(zind)];
vars3 = {'quint', 'wprime_quint',  'share_quint', 'renters_quint', 'a_quint', 'hp_quint_ho', 'ah_quint_ho', 'x_dhp_quint_ho', ...
    'hpminners_quint', 'high_risk_quint_ho', 'ah_quint_reg', 'x_dhp_quint_reg', 'F_ho_quint', 'x_dhp_quint_hpmin', ...
    'ah_quint_hpmin', 'a_quint', 'a_quint_hpmin', 'ah_quint_bho', 'hpminners_quint_ho', 'm_quint', 'hs_quint', ...
    'aid_quint_reg', 'subs_quint_reg', 'income_quint_reg', 'ah_hpmin', 'ah_ho', 'ah_bho', 'ah_quint_nhpmin', 'ah_nhpmin'};

exper = 'b';
eval(strcat('f_vec =f_vec_', exper, ';'));
eval(strcat('hp_vec =hp_vec_', exper, ';'));
eval(strcat('ha_vec =ha_vec_', exper, ';'));
eval(strcat('a_vec =a_vec_', exper, ';'));
eval(strcat('ah_vec =ah_vec_', exper, ';'));
eval(strcat('x_dhp_vec =x_dhp_vec_', exper, ';'));

wprime_vec = a_vec + hp_vec + ha_vec/(1+eta_b); %Next period's pre-shock wealth

%Assign indicies to quintiles
mat1 = [wprime_vec, f_vec];
[mat2, inds]  = sortrows(mat1, 1);
quint = zeros(length(f_vec), 1);
cum_sum =0;
for i = 1:length(f_vec)
    cum_sum = cum_sum + mat2(i, 2);
    if cum_sum <= 0.2 quint(inds(i)) =1;
    elseif cum_sum <= 0.4 quint(inds(i)) =2;
    elseif cum_sum <= 0.6 quint(inds(i)) = 3;
    elseif cum_sum <=0.8 quint(inds(i)) = 4;
    else quint(inds(i)) =5;
    end
end

%Vectors for aggregates
renters = zeros(size(hp_vec));
renters(hp_vec ==0) = 1;
hpminners = zeros(size(hp_vec));
hpminners(hp_vec == hpmin)=1;
high_risk_vec = zeros(size(hp_vec));
high_risk_vec(reg_vec ==2)=1;

%Averages by quintile
for i = 1:5
    wprime_quint(i) = f_vec(quint==i)'*wprime_vec(quint ==i)/sum(f_vec(quint==i));
    share_quint(i) = sum(f_vec(quint==i));
    renters_quint(i) = f_vec(quint==i)'*renters(quint ==i)/sum(f_vec(quint==i));
    a_quint(i) = f_vec(quint==i)'*a_vec(quint ==i)/sum(f_vec(quint==i));
    hpminners_quint(i) = f_vec(quint==i)'*hpminners(quint ==i)/sum(f_vec(quint==i));
    m_quint(i) = f_vec(quint==i)'*m_vec(quint ==i)/sum(f_vec(quint==i));
    hs_quint(i) = f_vec(quint==i)'*hs_vec(quint ==i)/sum(hs_vec(quint==i));
end


%Averages by quintile among homeowners
f_vec_ho = f_vec;
f_vec_ho(hp_vec ==0)=0;
for i = 1:5
    hp_quint_ho(i) = f_vec_ho(quint==i)'*hp_vec(quint ==i)/sum(f_vec_ho(quint==i));
    ah_quint_ho(i) = f_vec_ho(quint==i)'*ah_vec(quint ==i)/sum(f_vec_ho(quint==i));
    a_quint_ho(i) = f_vec_ho(quint==i)'*a_vec(quint ==i)/sum(f_vec_ho(quint==i));
    x_dhp_quint_ho(i) = f_vec_ho(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_ho(quint==i));
    hpminners_quint_ho(i) = f_vec_ho(quint==i)'*hpminners(quint ==i)/sum(f_vec_ho(quint==i));
    high_risk_quint_ho(i) = f_vec_ho(quint==i)'*high_risk_vec(quint ==i)/sum(f_vec_ho(quint==i));
    F_ho_quint(i) =  f_vec_ho(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_ho(quint==i));
end
ah_ho = f_vec_ho'*ah_vec/sum(f_vec_ho);

%Averages among homeowners at the mininium house size constraint
f_vec_hpmin = f_vec;
f_vec_hpmin(hp_vec > hpmin)=0;
f_vec_hpmin(hp_vec < hpmin) =0;
for i = 1:5
    hp_quint_hpmin(i) = f_vec_hpmin(quint==i)'*hp_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
    a_quint_hpmin(i) = f_vec_hpmin(quint==i)'*a_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
    ah_quint_hpmin(i) = f_vec_hpmin(quint==i)'*ah_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
    F_ho_quint_hpmin(i) =  f_vec_hpmin(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_hpmin(quint==i));
    x_dhp_quint_hpmin(i) = f_vec_hpmin(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
    hpminners_quint_hpmin(i) = f_vec_hpmin(quint==i)'*hpminners(quint ==i)/sum(f_vec_hpmin(quint==i));
    high_risk_quint_hpmin(i) = f_vec_hpmin(quint==i)'*high_risk_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
end
ah_hpmin = f_vec_hpmin'*ah_vec/sum(f_vec_hpmin);

%Averages among homeowners not the mininium house size constraint
f_vec_nhpmin = f_vec;
f_vec_nhpmin(hp_vec <= hpmin)=0;
for i = 1:5
    hp_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*hp_vec(quint ==i)/sum(f_vec_nhpmin(quint==i));
    a_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*a_vec(quint ==i)/sum(f_vec_nhpmin(quint==i));
    ah_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*ah_vec(quint ==i)/sum(f_vec_nhpmin(quint==i));
    F_ho_quint_nhpmin(i) =  f_vec_nhpmin(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_nhpmin(quint==i));
    x_dhp_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_nhpmin(quint==i));
    hpminners_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*hpminners(quint ==i)/sum(f_vec_nhpmin(quint==i));
    high_risk_quint_nhpmin(i) = f_vec_nhpmin(quint==i)'*high_risk_vec(quint ==i)/sum(f_vec_nhpmin(quint==i));
end
ah_nhpmin = f_vec_nhpmin'*ah_vec/sum(f_vec_nhpmin);


%Averages among homeowners barely above the minimum housesize
%constraint
f_vec_bho = f_vec;
f_vec_bho(hp_vec <= hpmin) =0; 
f_vec_bho(hp_vec >(1+x2)*hpmin)=0; 
for i = 1:5
    hp_quint_bho(i) = f_vec_bho(quint==i)'*hp_vec(quint ==i)/sum(f_vec_bho(quint==i));
    a_quint_bho(i) = f_vec_bho(quint==i)'*a_vec(quint ==i)/sum(f_vec_bho(quint==i));
    ah_quint_bho(i) = f_vec_bho(quint==i)'*ah_vec(quint ==i)/sum(f_vec_bho(quint==i));
    F_ho_quint_bho(i) =  f_vec_bho(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_bho(quint==i));
    x_dhp_quint_bho(i) = f_vec_bho(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_bho(quint==i));
end
ah_bho = f_vec_bho'*ah_vec/sum(f_vec_bho);


%Averages by quintile and region
for regi=1:2
    f_vec_ho_reg = f_vec_ho;
    f_vec_ho_reg(reg_vec ~=regi) =0;
    f_vec_reg = f_vec;
    f_vec_reg(reg_vec~=regi) =0;
    
    for i = 1:5
        ah_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*ah_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
        x_dhp_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
        hp_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*hp_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
        a_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*a_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
        aid_quint_reg(regi, i) = psiho_b*f_vec_ho_reg(quint==i)'*(Omegah_reg_b(regi)*exp(-theta*ah_vec(quint==i).^(1/theta)).*...
            hp_vec(quint ==i))/sum(f_vec_ho_reg(quint==i));
        subs_quint_reg(regi, i) =eta_b*f_vec_ho_reg(quint==i)'*ha_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
        income_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*w_reg_b(regi)*z_vec((quint==i))/sum(f_vec_ho_reg(quint==i));
    end
end

for i = 1:length(vars3)
    str = strcat(vars3{i}, '_', exper, '=', vars3{i}, ';');
    eval(str);
    str2 =  strcat(vars3{i}, '_', exper);
    if i ==1 save('quints', str2);
    else
        save('quints', str2, '-append');
    end
end

%% Aggregates by quintile in cc experiments

vars2 = {'f_vec', 'hp_vec', 'ha_vec', 'a_vec', 'ah_vec', 'x_dhp_vec', 'psiho', 'Omegah_reg', 'w_reg'};

for  k = 1:1:length(exper_vec2)
    exper = exper_vec2{k};
    for si = 0:1:nc
        for i = 1:length(vars2)
            eval(strcat(vars2{i}, '=', vars2{i}, '_', exper, num2str(si), ';'))
        end
        
        wprime_vec = a_vec + hp_vec + ha_vec/(1+eta_b); %Next period's pre-shock wealth
        wprime_vec = a_vec + hp_vec + ha_vec/(1+eta_b); %Next period's pre-shock wealth
        
        %Assign indicies to quintiles
        mat1 = [wprime_vec, f_vec];
        [mat2, inds]  = sortrows(mat1, 1);
        quint = zeros(length(f_vec), 1);
        cum_sum =0;
        for i = 1:length(f_vec)
            cum_sum = cum_sum + mat2(i, 2);
            if cum_sum <= 0.2 quint(inds(i)) =1;
            elseif cum_sum <= 0.4 quint(inds(i)) =2;
            elseif cum_sum <= 0.6 quint(inds(i)) = 3;
            elseif cum_sum <=0.8 quint(inds(i)) = 4;
            else quint(inds(i)) =5;
            end
        end
        
        %Vectors for aggregates
        renters = zeros(size(hp_vec));
        renters(hp_vec ==0) = 1;
        hpminners = zeros(size(hp_vec));
        hpminners(hp_vec == hpmin)=1;
        high_risk_vec = zeros(size(hp_vec));
        high_risk_vec(reg_vec ==2)=1;
        
        %Averages by quintile
        for i = 1:5
            wprime_quint(i) = f_vec(quint==i)'*wprime_vec(quint ==i)/sum(f_vec(quint==i));
            share_quint(i) = sum(f_vec(quint==i));
            renters_quint(i) = f_vec(quint==i)'*renters(quint ==i)/sum(f_vec(quint==i));
            a_quint(i) = f_vec(quint==i)'*a_vec(quint ==i)/sum(f_vec(quint==i));
            hpminners_quint(i) = f_vec(quint==i)'*hpminners(quint ==i)/sum(f_vec(quint==i));
            m_quint(i) = f_vec(quint==i)'*m_vec(quint ==i)/sum(f_vec(quint==i));
            hs_quint(i) = f_vec(quint==i)'*hs_vec(quint ==i)/sum(hs_vec(quint==i));
            num_quint(i) = sum(f_vec(quint==i));
        end
        
        %Averages by quintile among homeowners
        f_vec_ho = f_vec;
        f_vec_ho(hp_vec < hpmin) =0;
        for i = 1:5
            hp_quint_ho(i) = f_vec_ho(quint==i)'*hp_vec(quint ==i)/sum(f_vec_ho(quint==i));
            a_quint_ho(i) = f_vec_ho(quint==i)'*a_vec(quint ==i)/sum(f_vec_ho(quint==i));
            ah_quint_ho(i) = f_vec_ho(quint==i)'*ah_vec(quint ==i)/sum(f_vec_ho(quint==i));
            F_ho_quint_ho(i) =  f_vec_ho(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_ho(quint==i));
            x_dhp_quint_ho(i) = f_vec_ho(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_ho(quint==i));
            hpminners_quint_ho(i) = f_vec_ho(quint==i)'*hpminners(quint ==i)/sum(f_vec_ho(quint==i));
            high_risk_quint_ho(i) = f_vec_ho(quint==i)'*high_risk_vec(quint ==i)/sum(f_vec_ho(quint==i));
        end
        ah_ho = f_vec_ho'*ah_vec/sum(f_vec_ho);
        
        
        %Averages among homeowners at the mininium house size constraint
        f_vec_hpmin = f_vec;
        f_vec_hpmin(hp_vec > hpmin)=0;
        f_vec_hpmin(hp_vec < hpmin) =0;
        for i = 1:5
            hp_quint_hpmin(i) = f_vec_hpmin(quint==i)'*hp_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
            a_quint_hpmin(i) = f_vec_hpmin(quint==i)'*a_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
            ah_quint_hpmin(i) = f_vec_hpmin(quint==i)'*ah_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
            F_ho_quint_hpmin(i) =  f_vec_hpmin(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_hpmin(quint==i));
            x_dhp_quint_hpmin(i) = f_vec_hpmin(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
            hpminners_quint_hpmin(i) = f_vec_hpmin(quint==i)'*hpminners(quint ==i)/sum(f_vec_hpmin(quint==i));
            high_risk_quint_hpmin(i) = f_vec_hpmin(quint==i)'*high_risk_vec(quint ==i)/sum(f_vec_hpmin(quint==i));
        end
        ah_hpmin = f_vec_hpmin'*ah_vec/sum(f_vec_hpmin);
        
        %Averages among homeowners barely above the minimum housesize
        %constraint
        f_vec_bho = f_vec;
        f_vec_bho(hp_vec <= hpmin) =0;
        f_vec_bho(hp_vec >(1+x2)*hpmin)=0;
        for i = 1:5
            hp_quint_bho(i) = f_vec_bho(quint==i)'*hp_vec(quint ==i)/sum(f_vec_bho(quint==i));
            a_quint_bho(i) = f_vec_bho(quint==i)'*a_vec(quint ==i)/sum(f_vec_bho(quint==i));
            ah_quint_bho(i) = f_vec_bho(quint==i)'*ah_vec(quint ==i)/sum(f_vec_bho(quint==i));
            F_ho_quint_bho(i) =  f_vec_bho(quint==i)'*exp(-theta*ah_vec(quint ==i).^(1/theta))/sum(f_vec_bho(quint==i));
            x_dhp_quint_bho(i) = f_vec_bho(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_bho(quint==i));
        end
        ah_bho = f_vec_bho'*ah_vec/sum(f_vec_bho);
        
        
        %Averages by quintile and region
        for regi=1:2
            f_vec_ho_reg = f_vec_ho;
            f_vec_ho_reg(reg_vec ~=regi) =0;
            f_vec_reg = f_vec;
            f_vec_reg(reg_vec~=regi) =0;
            
            for i = 1:5
                ah_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*ah_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
                x_dhp_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*x_dhp_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
                hp_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*hp_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
                a_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*a_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
                renters_quint_reg(regi, i) = f_vec_reg(quint==i)'*renters(quint ==i)/sum(f_vec_reg(quint==i));
                aid_quint_reg(regi, i) = psiho*f_vec_ho_reg(quint==i)'*(Omegah_reg(regi)*exp(-theta*ah_vec(quint==i).^(1/theta)).*...
                    hp_vec(quint ==i))/sum(f_vec_ho_reg(quint==i));
                subs_quint_reg(regi, i) =eta*f_vec_ho_reg(quint==i)'*ha_vec(quint ==i)/sum(f_vec_ho_reg(quint==i));
                income_quint_reg(regi, i) = f_vec_ho_reg(quint==i)'*w_reg(regi)*z_vec((quint==i))/sum(f_vec_ho_reg(quint==i));
            end
        end
        
        for i = 1:length(vars3)
            str = strcat(vars3{i}, '_', exper, num2str(si), '=', vars3{i}, ';');
            eval(str);
            str2 =  strcat(vars3{i}, '_', exper, num2str(si));
            save('quints', str2, '-append');
        end
    end
end




