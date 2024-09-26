%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: tables.m
% Bay: Stephie Fried
% Date: Spring 2018
% Purpose: Calculates tables and figures for model results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all
bpath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output';
tablepath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output\tables';
figpath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output\figures';
cd(bpath);

%% Parameter Tables
load params_targs_base
load agg_b

cd(tablepath);
FID = fopen('params_direct.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Parameter Values: External Calibration} \\label{tab:params_direct}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l l}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Parameter & Value  \\\\ \n ');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Low-risk region storm probability: $\\gamma_1$ &  %8.3f   \\\\ ', rho_l);
fprintf(FID, ' High-risk region storm probability: $\\gamma_2$ &  %8.3f  \\\\ ', rho_h);
fprintf(FID, ' Coefficient of relative risk aversion: $\\sigma$ &  %8.0f    \\\\ ', sigma);
fprintf(FID, ' Depreciation rate of final-good capital: $\\delta^y$ &  %8.2f    \\\\ ', deltak);
fprintf(FID, ' Depreciation rate of housing capital: $\\delta^h$ &  %8.2f   \\\\ ', deltah);
fprintf(FID, ' Depreciation rate of adaptation capital: $\\delta^a$ &  %8.2f    \\\\ ', deltaa);
fprintf(FID, ' TFP in final-good production: $A^y$ &  %8.0f  \\\\ ', 1);
fprintf(FID, ' TFP in housing services: $A^h$ &  %8.2f \\\\ ', Ah);
fprintf(FID, ' Capital share in final-good production: $\\alpha$ &  %8.2f   \\\\ ', alpha);
fprintf(FID, 'International interest rate: $r^\\star$ &  %8.2f \\\\ ', r);
fprintf(FID, 'Persistent shock persistence: $\\rho$ &  %8.2f  \\\\ ', 0.97);
fprintf(FID, 'Persistent shock innovation variance: $\\sigma_{\\eps}^2$ &  %8.2f   \\\\ ', 0.016);
fprintf(FID, 'Fixed-effect variance: $\\sigma_{\\xi}^2$ &  %8.2f    \\\\ ', 0.66);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

cd(tablepath);
FID = fopen('params_mom.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Parameter Values: Internal Calibration} \\label{tab:params_mom}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l l}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Parameter & Value  \\\\ \n ');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Low-risk, housing, storm severity: $\\Omega_1^h$ &  %8.2f \\\\ ', Omegah_reg_b(1));
fprintf(FID, ' High-risk, housing, storm severity: $\\Omega_2^h$ &  %8.2f \\\\ ', Omegah_reg_b(2));
fprintf(FID, ' Low-risk, final-good, storm severity: $\\Omega_1^y$ &  %8.2f \\\\ ', Omegak_reg_b(1));
fprintf(FID, ' High-risk, final-good, storm severity: $\\Omega_2^y$ &  %8.2f \\\\ ', Omegak_reg_b(2));
fprintf(FID, ' Effectiveness of adaptation: $\\theta$ &  %8.2f   \\\\ ', theta);
fprintf(FID, ' Subsidy for adaptation investment: $\\eta$ &  %8.2f  \\\\ ', eta);
fprintf(FID, ' Disaster aid fraction for homeowners: $\\psi^{ho}$ &  %8.2f   \\\\ ', psiho);
fprintf(FID, ' Disaster aid fraction for renters: $\\psi^{hr}$ &  %8.2f  \\\\ ', psihr);
fprintf(FID, ' Disaster aid fraction for rental-housing firms: $\\psi^{kr}$ &  %8.2f  \\\\ ', psikr);
fprintf(FID, ' Disaster aid fraction  for final-good firms: $\\psi^{ky}$ &  %8.2f  \\\\ ', psiky);
fprintf(FID, ' Discount factor: $\\beta$ &  %8.2f  \\\\ ', beta);
fprintf(FID, ' Non-durable consumption exponent: $\\zeta$ &  %8.2f  \\\\ ', zeta);
fprintf(FID, ' Preference discount from renting: $\\bar{\\phi}^r$ &  %8.2f   \\\\ ', phi_r);
fprintf(FID, ' Preference discount from owning: $\\bar{\\phi}^o$ &  <1  \\\\ ');
fprintf(FID, ' Fraction of HHs that prefer homeownership: $\\chi$ &  %8.2f  \\\\ ', chi);
fprintf(FID, ' Disbursement cost of insurance: $\\lambda$ &  %8.3f  \\\\ ', lambda);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% %Model fit tables

ffaid_aid_total_aid = 1 - renter_aid_total_aid - owner_aid_total_aid - rental_housing_aid_total_aid; 
ffaid_aid_total_aid_target = 1 - renter_aid_total_aid_target - owner_aid_total_aid_target - rental_housing_aid_total_aid_target; 

%Targeted moments
cd(tablepath);
FID = fopen('targets.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Targeted Moments} \\label{tab:targets}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' Moment & Model value & Empirical value  \\\\ \n ');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Relative aid  &  %8.2f  &  %8.2f \\\\ ', rel_aid, rel_aid_target);
fprintf(FID, 'Relative severity  &  %8.2f  & %8.2f  \\\\ ', Omegah_reg_b(2)/Omegah_reg_b(1), Omegah_reg_b(2)/Omegah_reg_b(1));
fprintf(FID, '(Disaster aid)/GDP &  %8.1e  &  %8.1e \\\\ ', aid_output, aid_output_target);
fprintf(FID, '(Adaptation subsidy)/GDP &  %8.1e  &  %8.1e  \\\\ ', subs_output, subs_output_target);
fprintf(FID, '(Aid to homeowners)/(total disaster aid) &  %8.2f  &   %8.2f \\\\ ', owner_aid_total_aid, owner_aid_total_aid_target);
fprintf(FID, '(Aid to renters)/(total disaster aid) &  %8.2f  &  %8.2f \\\\ ', renter_aid_total_aid, renter_aid_total_aid_target);
fprintf(FID, '(Aid to final-good firms)/(total disaster aid) &  %8.2f  &  %8.2f \\\\ '...
, ffaid_aid_total_aid, ffaid_aid_total_aid_target);
fprintf(FID, '(Disaster aid)/(storm damage) &  %8.2f  &  %8.2f \\\\ '...
, fema_damage, fema_damage_target);
fprintf(FID, '(Damaged housing capital)/(damaged capital) &  %8.2f  &  %8.2f \\\\ '...
, (Dhp_b + Dhr_b)/D_cap_b, 0.43);
fprintf(FID, '(Insurance claims)/(owner-occupied housing damage) &  %8.2f  & %8.2f \\\\ ', insurance_house_damage, insurance_house_damage_target);
fprintf(FID, '(Residential assets)/(non-residential assets) &  %8.2f  &  %8.2f \\\\ ', housing_capital, housing_capital_target);
fprintf(FID, '(Owner-occupied housing)/(total housing) &  %8.2f  &  %8.2f \\\\ ', owned_total, owned_total_target);
fprintf(FID, '(Avg. income homeowners)/(avg. income renters) &  %8.1f  &  %8.1f  \\\\ ', owner_renter_inc, owner_renter_inc_target);
fprintf(FID, 'Homeowner fraction &  %8.2f  & %8.2f \\\\ ', owner_frac, owner_frac_target);
fprintf(FID, '(Net wealth)/GDP  &  %8.1f  &  %8.1f  \\\\ ', wealth_output, wealth_output_target);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Income and wealth shares by quintile

%Income quintile table
cd(tablepath);
FID = fopen('income_quintile.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Share of Income Received by Each Quintile of the Income Distribution} \\label{tab:income_quintile}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '& \\multicolumn{5}{c}{Quintile} \\\\ \n');
fprintf(FID, '\\cline{2-6}  \n');
fprintf(FID, '  & First & Second & Third & Fourth & Fifth   \\\\ \n ');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Model &  %8.1f  &  %8.1f &  %8.1f  &  %8.1f &  %8.1f   \\\\ ',i1, i2, i3, i4, i5);
fprintf(FID, 'Data &  %8.1f  &  %8.1f &  %8.1f  &  %8.1f &  %8.1f   \\\\ ',i1_target, i2_target, i3_target, i4_target, i5_target);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%Wealth quintile table
cd(tablepath);
FID = fopen('wealth_quintile.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Share of Wealth Held by Each Quintile of the Wealth Distribution} \\label{tab:wealth_quintile}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '& \\multicolumn{5}{c}{Quintile} \\\\ \n');
fprintf(FID, '\\cline{2-6}  \n');
fprintf(FID, '  & First & Second & Third & Fourth & Fifth   \\\\ \n ');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Model &  %8.1f  &  %8.1f &  %8.1f  &  %8.1f &  %8.1f   \\\\ ',w1, w2, w3, w4, w5);
fprintf(FID, 'Data &  %8.1f  &  %8.1f &  %8.1f  &  %8.1f &  %8.1f   \\\\ ',w1_target, w2_target, w3_target, w4_target, w5_target);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Tables: 5.1

%Baseline aggregates
load agg_b
cd(tablepath);
FID = fopen('baseline_agg.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Macro-Aggregates: Baseline} \\label{tab:baseline_agg}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& Low-risk & High-risk & \\multirow{2}{*}{Aggregate} \\\\[-0.75ex] \n');
fprintf(FID, ' 	& region & region &  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Adaptation capital & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f}  \\\\[-0.75ex]', ...
    frac_adapt_reg_b(1)*100, frac_adapt_reg_b(2)*100, frac_adapt_b*100 );
fprintf(FID, '(percent of total capital) \\\\');
fprintf(FID, 'Damaged capital & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f}  \\\\[-0.75ex]', ...
   D_totalK_reg_b(1)/total_cap_reg_b(1)*100,D_totalK_reg_b(2)/total_cap_reg_b(2)*100, D_cap_b/total_cap_b*100 );
fprintf(FID, '(percent of total capital) \\\\');
fprintf(FID, 'Damage & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f}  \\\\[-0.75ex]', ...
   D_reg_b(1)/Ytilde_reg_b(1)*100,D_reg_b(2)/Ytilde_reg_b(2)*100, D_b/Ytilde_b*100 );
fprintf(FID, '(percent of output) \\\\');
fprintf(FID, 'Insurance claims & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f}  \\\\[-0.75ex]', ...
   x_dhp_reg_b(1)*100, x_dhp_reg_b(2)*100, x_dhp_b*100 );
fprintf(FID, '(percent of damaged owner-occupied housing capital) \\\\');
fprintf(FID, 'Disaster aid + subsidy & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f}  \\\\[-0.75ex]', ...
   (aid_reg_b(1) + subs_reg_b(1))/(tau_b*w_reg_b(1))*100,(aid_reg_b(2) + subs_reg_b(2))/(tau_b*w_reg_b(2))*100, 100.0 );
fprintf(FID, '(percent of model tax payments) \\\\');
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%Adaptive capacity
load agg_b
cd(tablepath);
FID = fopen('adaptive_capacity.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Adaptive Capacity} \\label{tab:adaptive_capacity}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& Low-risk & High-risk & \\multirow{2}{*}{Aggregate} \\\\[-0.75ex] \n');
fprintf(FID, ' 	& region & region &  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Homeowners: $h^a/h^p$ & %8.1s & %8.1s & %8.1s  \\\\[0.5ex]', ...
   ah_reg_b(1), ah_reg_b(2), ah_b );
fprintf(FID, 'Rental-housing firms: $k^{ra}/k^{rp}$ & %8.1s & %8.1s & %8.1s  \\\\[0.5ex]', ...
   ar_reg_b(1), ar_reg_b(2), mean(ar_reg_b) );
fprintf(FID, 'Final-good firms: $k^{ya}/k^{yp}$ & %8.1s & %8.1s & %8.1s  \\\\[0.5ex]', ...
   ay_reg_b(1), ay_reg_b(2), mean(ay_reg_b) );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%Adaptation and insurance by quintile
load quints
cd(tablepath);
FID = fopen('quintiles.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Homeowner Adaptive Capacity and Insurance By Wealth Quintile} \\label{tab:quintiles}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' & \\multicolumn{5}{c}{Quintile}\\\\');
fprintf(FID, '  \\cline{2-6}');
fprintf(FID, ' 	&  1 & 2 &  3&  4 &  5 \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Insurance claims & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f}  & \\multirow{2}{*}{%8.1f} & \\multirow{2}{*}{%8.1f} \\\\[-0.75ex]', ...
   x_dhp_quint_ho_b*100'  );
fprintf(FID, '(percent of damaged housing capital) \\\\');
fprintf(FID, 'Adaptive capacity: all homeowners & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f}  \\\\', ...
   ah_quint_ho_b'  );
fprintf(FID, 'Adaptive capacity: $h^p = \\underline{h}^p$ & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{%8.3f} & \\multirow{1}{*}{-}  \\\\', ...
   ah_quint_hpmin_b(1:4)'  );
%fprintf(FID, 'Percent of homeowners with $h^p = \\underline{h}^p$ & \\multirow{1}{*}{%8.1f} & \\multirow{1}{*}{%8.1f} & \\multirow{1}{*}{%8.1f}  & \\multirow{1}{*}{%8.1f} & \\multirow{1}{*}{%8.1f} \\\\', ...
 %  hpminners_quint*100'  );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Tables 5.2: Federal policy experiments
load agg_b
load agg_np
load agg_ns

%Effects of federal policy on adaptatation 
cd(tablepath);
FID = fopen('federal_policy.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Effects of Federal Policy on Adaptation, Storm Damage, and Insurance} \\label{tab:federal_policy}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& No policy & Disaster aid only & Baseline \\\\[-0.5ex] \n');
fprintf(FID, ' 	& ($\\psi = \\eta =0$) & ($\\psi >0, \\; \\eta =0$)  &  ($\\psi >0, \\; \\eta >0$)  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Adaptation capital &  \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f}  \\\\[-0.75ex]',...
    frac_adapt_np*100, frac_adapt_ns*100, frac_adapt_b*100);
fprintf(FID, '(percent of total capital) \\\\');
fprintf(FID, 'Damage & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f} & \\multirow{2}{*}{%8.2f}  \\\\[-0.75ex]',...
    (D_np/Ytilde_np)*100, (D_ns/Ytilde_ns)*100, (D_b/Ytilde_b)*100);
fprintf(FID, '(percent of output) \\\\');
fprintf(FID, 'Insurance claims & \\multirow{3}{*}{%8.1f} & \\multirow{3}{*}{%8.1f} & \\multirow{3}{*}{%8.1f}  \\\\[-0.75ex]', ...
   x_dhp_np*100,x_dhp_ns*100, x_dhp_b*100 );
fprintf(FID, '(percent of damaged owner- \\\\[-0.75ex]');
fprintf(FID, 'occupied housing capital) \\\\');
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

load value_policy_funcs
%Welfare effects of federal policy
cd(tablepath);
FID = fopen('fed_welfare.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Welfare Effects of Federal Disaster Policy: CEV (percent)} \\label{tab:fed_welfare}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& Low-risk & High-risk & \\multirow{2}{*}{Aggregate} \\\\[-0.75ex] \n');
fprintf(FID, ' 	& region & region &  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Aid only ($\\psi >0, \\; \\eta =0$) &  \\multirow{1}{*}{%8.2f} & \\multirow{1}{*}{%8.2f} & \\multirow{1}{*}{%8.3f}  \\\\'...
    , ((wel_l_ns/wel_l_np)^(1/(zeta*(1-sigma)))-1)*100, ((wel_h_ns/wel_h_np)^(1/(zeta*(1-sigma)))-1)*100,...
    ((wel_ns/wel_np)^(1/(zeta*(1-sigma)))-1)*100);
fprintf(FID, 'Baseline ($\\psi >0, \\; \\eta >0$) &  \\multirow{1}{*}{%8.2f} & \\multirow{1}{*}{%8.2f} & \\multirow{1}{*}{%8.3f}  \\\\'...
    , ((wel_l_b/wel_l_np)^(1/(zeta*(1-sigma)))-1)*100, ((wel_h_b/wel_h_np)^(1/(zeta*(1-sigma)))-1)*100,...
    ((wel_b/wel_np)^(1/(zeta*(1-sigma)))-1)*100);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%% Climate change experiments
load agg_b
load agg_cc0
load agg_cc1
load agg_cc2
load agg_cc3
load agg_cc4
load agg_cc5
load agg_ccna0
load agg_ccna1
load agg_ccna2
load agg_ccna3
load agg_ccna4
load agg_ccna5
load agg_ccna9

%Parameter value table
cd(tablepath);
FID = fopen('climate_params.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Parameter Values for Climate Change Experiments} \\label{tab:climate_params}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '& \\multicolumn{2}{c}{Storm probability} & \\multicolumn{4}{c}{Storm Severity} \\\\ \n');
fprintf(FID, '		\\cline{2-3} \\cline{4-7} \n');
fprintf(FID, '		&&&&&& \\\\[-2.5ex]  \n');
fprintf(FID, ' 	& $\\gamma_1$ & $\\gamma_2$&  $\\Omega^h_1$ & $\\Omega^y_1$ & $\\Omega^h_2$ & $\\Omega^y_2$ \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Baseline & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_b(1),rho_reg_b(2), Omegah_reg_b(1), ...
    Omegak_reg_b(1),Omegah_reg_b(2), Omegak_reg_b(2) );
fprintf(FID, 'Cyclone method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc0(1),rho_reg_cc0(2), ...
    Omegah_reg_cc0(1), Omegak_reg_cc0(1),Omegah_reg_cc0(2), Omegak_reg_cc0(2) );
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc1(1),rho_reg_cc1(2), ...
    Omegah_reg_cc1(1), Omegak_reg_cc1(1),Omegah_reg_cc1(2), Omegak_reg_cc1(2) );
fprintf(FID, 'Extreme-precipitation method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc2(1),rho_reg_cc2(2), ...
    Omegah_reg_cc2(1), Omegak_reg_cc2(1),Omegah_reg_cc2(2), Omegak_reg_cc2(2) );
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc3(1),rho_reg_cc3(2), ...
    Omegah_reg_cc3(1), Omegak_reg_cc3(1),Omegah_reg_cc3(2), Omegak_reg_cc3(2) );
fprintf(FID, 'Flood method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc4(1),rho_reg_cc4(2), ...
    Omegah_reg_cc4(1), Omegak_reg_cc4(1),Omegah_reg_cc4(2), Omegak_reg_cc4(2) );
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.3f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\', rho_reg_cc5(1),rho_reg_cc5(2), ...
    Omegah_reg_cc5(1), Omegak_reg_cc5(1),Omegah_reg_cc5(2), Omegak_reg_cc5(2) );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);

%Adaptation capital table
cd(tablepath);
FID = fopen('adaptation_capital.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Adaptation Capital as a Percent of the Total Capital Stock} \\label{tab:adaptation_capital}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& Low-risk & High-risk & \\multirow{2}{*}{Aggregate} \\\\[-0.5ex] \n');
fprintf(FID, ' 	& region & region &  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Baseline & %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_b(1)*100, frac_adapt_reg_b(2)*100, frac_adapt_b*100 );
fprintf(FID, 'Cyclone method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5  & %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc0(1)*100, ...
    frac_adapt_reg_cc0(2)*100, frac_adapt_cc0*100 );
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc1(1)*100, ...
    frac_adapt_reg_cc1(2)*100, frac_adapt_cc1*100 );
fprintf(FID, 'Extreme-precipitation method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 & %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc2(1)*100, ...
    frac_adapt_reg_cc2(2)*100, frac_adapt_cc2*100 );
fprintf(FID, '\\hspace{1em} RCP 8.5& %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc3(1)*100, ...
    frac_adapt_reg_cc3(2)*100, frac_adapt_cc3*100 );
fprintf(FID, 'Flood method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5& %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc4(1)*100, ...
    frac_adapt_reg_cc4(2)*100, frac_adapt_cc4*100 );
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f & %8.2f  \\\\', frac_adapt_reg_cc5(1)*100, ...
    frac_adapt_reg_cc5(2)*100, frac_adapt_cc5*100 );
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


load value_policy_funcs
%Welfare effects of climate change
cd(tablepath);
FID = fopen('welfare_cc.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Welfare Effects of Climate Change With and Without Adaptation: CEV (percent)} \\label{tab:welfare_cc}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	& With adaptation & Without adaptation  \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Cyclone method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5  & %8.2f & %8.2f   \\\\', ((wel_cc0/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna0/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f   \\\\', ((wel_cc1/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna1/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, 'Extreme-precipitation method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 &  %8.2f & %8.2f   \\\\', ((wel_cc2/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna2/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5&  %8.2f & %8.2f   \\\\', ((wel_cc3/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna3/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, 'Flood method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5& %8.2f & %8.2f   \\\\', ((wel_cc4/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna4/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f   \\\\', ((wel_cc5/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccna5/wel_ccna9)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


val_b = [renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_b, ah_hpmin_b, ah_bho_b];

val_cc0 = [renters_reg_cc0/sum(renters_reg_cc0)*ar_reg_cc0'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc0-ah_b, ...
    ah_hpmin_cc0 - ah_hpmin_b, ah_bho_cc0- ah_bho_b];
val_cc1 = [renters_reg_cc1/sum(renters_reg_cc1)*ar_reg_cc1'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc1-ah_b,  ...
    ah_hpmin_cc1 - ah_hpmin_b, ah_bho_cc1- ah_bho_b];
val_cc2 = [renters_reg_cc2/sum(renters_reg_cc2)*ar_reg_cc2'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc2-ah_b,  ...
    ah_hpmin_cc2 - ah_hpmin_b, ah_bho_cc2- ah_bho_b];
val_cc3 = [renters_reg_cc3/sum(renters_reg_cc3)*ar_reg_cc3'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc3-ah_b,  ...
    ah_hpmin_cc3 - ah_hpmin_b, ah_bho_cc3- ah_bho_b];
val_cc4 = [renters_reg_cc4/sum(renters_reg_cc4)*ar_reg_cc4'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc4-ah_b, ...
    ah_hpmin_cc4 - ah_hpmin_b, ah_bho_cc4- ah_bho_b];
val_cc5 = [renters_reg_cc5/sum(renters_reg_cc5)*ar_reg_cc5'-renters_reg_b/sum(renters_reg_b)*ar_reg_b', ah_cc5-ah_b, ...
    ah_hpmin_cc5 - ah_hpmin_b, ah_bho_cc5- ah_bho_b];

cd(tablepath);
FID = fopen('adapt_cc.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Increase in Adaptive Capacity From the Baseline} \\label{tab:adapt_cc}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' &  & \\multicolumn{3}{c}{Homeowners}\\\\');
fprintf(FID, '  \\cline{3-5}');
fprintf(FID, ' 	& Renters &  \\multirow{2}{*}{All} & \\multirow{2}{*}{Constrained} &  Barely  \\\\[-1.25ex] \n');
fprintf(FID, ' 	& &  & &    unconstrained  \\\\ \n');
fprintf(FID, '\\hline \n');
%fprintf(FID, '\\Baseline  & %8.3f & %8.3f  & %8.3f & %8.3f     \\\\', val_cc0);
fprintf(FID, 'Cyclone method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5  & %8.3f & %8.3f  & %8.3f & %8.3f     \\\\', val_cc0);
fprintf(FID, '\\hspace{1em} RCP 8.5  & %8.3f & %8.3f  & %8.3f & %8.3f      \\\\', val_cc1);
fprintf(FID, 'Extreme-precipitation method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5   & %8.3f & %8.3f  & %8.3f & %8.3f   \\\\',val_cc2);
fprintf(FID, '\\hspace{1em} RCP 8.5  & %8.3f & %8.3f  & %8.3f & %8.3f  \\\\', val_cc3);
fprintf(FID, 'Flood method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5  & %8.3f & %8.3f  & %8.3f & %8.3f    \\\\', val_cc4);
fprintf(FID, '\\hspace{1em} RCP 8.5   & %8.3f & %8.3f  & %8.3f & %8.3f    \\\\', val_cc5);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);


%% Effect of adaptation on damage from climate change
pic = 1; 
somenames2={'Cylcone'; 'Extreme-precipitation'; 'Flood'};

%Percent increase in damage under RCP 4.5
plotd_45 = ([D_ccna0/D_ccna9, D_ccna2/D_ccna9, D_ccna4/D_ccna9; D_cc0/D_b, D_cc2/D_b, D_cc4/D_b]'-1)*100;

figure(1)
b = bar(plotd_45, 'Edgecolor', 'flat');
b(1).FaceColor=[0 0 0] + 0.75;
b(1).EdgeColor=[0 0 0] + 0.75;
b(2).FaceColor=[0 0 0] + 0.25;
b(2).EdgeColor=[0 0 0]+0.25;
xt = get(gca, 'XTick');
labels = arrayfun(@(value) num2str(value,'%2.1f'),plotd_45(:,1),'UniformOutput',false);
text(xt-.15,plotd_45(:,1), labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
labels = arrayfun(@(value) num2str(value,'%2.1f'),plotd_45(:,2),'UniformOutput',false);
text(xt+.15,plotd_45(:,2), labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('Method of Modeling Climate Change')
ylabel('Percent Increase in Damage')
set(gca, 'xticklabel', somenames2)
set(gca,'FontSize',12)
title('RCP 4.5')
ylim([0, 150])
legend('No Adaptation','Adaptation' ,'Location', 'NW')
legend boxoff
box off
if pic ==1
    saveas(gcf,[figpath ,filesep, 'D_45'],'eps');
end

%Percent increase in damage under RCP 8.5
plotd_85 = ([D_ccna1/D_b, D_ccna3/D_b, D_ccna5/D_b; D_cc1/D_b, D_cc3/D_b, D_cc5/D_b]'-1)*100;

figure(2)
b =bar(plotd_85, 'Edgecolor', 'flat');
b(1).FaceColor=[0 0 0] + 0.75;
b(1).EdgeColor=[0 0 0] + 0.75;
b(2).FaceColor=[0 0 0] + 0.25;
b(2).EdgeColor=[0 0 0]+0.25;
xt = get(gca, 'XTick');
labels = arrayfun(@(value) num2str(value,'%2.1f'),plotd_85(:,1),'UniformOutput',false);
text(xt-.15,plotd_85(:,1), labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
labels = arrayfun(@(value) num2str(value,'%2.1f'),plotd_85(:,2),'UniformOutput',false);
text(xt+.15,plotd_85(:,2), labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('Method of Modeling Climate Change')
ylabel('Percent Increase in Damage')
set(gca, 'xticklabel', somenames2)
set(gca,'FontSize',12)
title('RCP 8.5')
%ylim([0, 70])
legend('No Adaptation','Adaptation' ,'Location', 'NW')
legend boxoff
box off
if pic ==1
    saveas(gcf,[figpath ,filesep, 'D_85'],'eps');
end

%% 5.4 Idiosyncratic Risk

load value_policy_funcs
%Welfare effects of climate change
cd(tablepath);
FID = fopen('welfare_ccnr.tex', 'w');
fprintf(FID, '\\begin{table}[H] \n');
fprintf(FID, '\\singlespace \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\caption{Welfare Effects of Climate Change: \\\\ Idiosyncratic Storm Risk Compared to Damage-Function Specifications} \\label{tab:welfare_ccnr}\n');
fprintf(FID, '\\vspace{-.1in} \n');
fprintf(FID, '\\begin{tabular}{l c c}\\hline \n');
fprintf(FID, '\\hline \n');
fprintf(FID, ' 	&Idiosyncratic Storm Risk & Damage Function \\\\ \n');
fprintf(FID, '\\hline \n');
fprintf(FID, 'Cyclone method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5  & %8.2f & %8.2f   \\\\', ((wel_cc0/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr0/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f   \\\\', ((wel_cc1/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr1/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, 'Extreme-precipitation method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5 &  %8.2f & %8.2f   \\\\', ((wel_cc2/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr2/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5&  %8.2f & %8.2f   \\\\', ((wel_cc3/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr3/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, 'Flood method \\\\ \n');
fprintf(FID, '\\hspace{1em} RCP 4.5& %8.2f & %8.2f   \\\\', ((wel_cc4/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr4/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hspace{1em} RCP 8.5 & %8.2f & %8.2f   \\\\', ((wel_cc5/wel_b)^(1/(zeta*(1-sigma))) -1)*100, ...
    ((wel_ccnr5/wel_nr)^(1/(zeta*(1-sigma))) -1)*100);
fprintf(FID, '\\hline \n');
fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{table} \n');
fclose(FID);
cd(bpath);









  
  