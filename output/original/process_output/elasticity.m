%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: elasticity.m
% Bay: Stephie Fried
% Date: Summer 2021 
% Purpose: Calculates disaster aid and damage elasticities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all

datapath = strcat('C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\c_programs\results_1109');
bpath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output';
cd(bpath);

load agg_hr
load agg_b

%Damage elasticity
top1 = (Dk_reg_hr + Dhp_reg_hr + Dhr_reg_hr)./(rho_reg_hr.*(GDP_reg_hr));
top2 =  (Dk_reg_b + Dhp_reg_b + Dhr_reg_b)./(rho_reg_b.*(GDP_reg_b));
top = top1./top2 -1;
bot =(rho_reg_hr - rho_reg_b)./rho_reg_b; 

damage_elas_reg = top./bot *100;
damage_elas = mean(damage_elas_reg)

%Aid elasticity
top1 = (aid_reg_hr + aid_reg_hr + aid_reg_hr)./(rho_reg_hr.*(GDP_reg_hr));
top2 =  (aid_reg_b + aid_reg_b + aid_reg_b)./(rho_reg_b.*(GDP_reg_b));
top = top1./top2 -1;
bot =(rho_reg_hr - rho_reg_b)./rho_reg_b; 

aid_elas_reg = top./bot *100;
aid_elas = mean(aid_elas_reg)