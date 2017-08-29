library(deSolve)
library(reshape2)
####################### PBPK model #################################

# Inhalation exposure during labor
# Experimental setup: 2h work followed by  5 hour follow-up
# Source: Ernsgaard, 2014
# R script: by Signe Klinting, 2017

####################### PBPK model #################################

pbpkEx <- function(Time, Concentration, parameters){
  with(as.list(c(Concentration, parameters)), {
    AW 		= ifelse(Time %% (24*60) >= (2*60), 0, 1)	
    Cair 		= ppm/molvol*AW				 		
    Qp 		= ifelse(AW == 1, Qp_work, Qp_rest)		
    W 		= AW*50								
    
    VOD 		= 10.1*W/1000								
    Qfat 		= 0.022*Vfat + VOD						
    Qkidney 	= 0.0452*Vkidney*(100-14*VOD)				
    Qvrg 		= Qbrain + Qkidney + Qother				
    Qwm 		= Qrm + 7.38*VOD						
    Qliv 		= (0.00271+0.00833)*Vliv*(100-19*VOD)		
    Qalv 		= QaQp*Qp								
    Qtot 		= Qrm + Qwm + Qfat + Qliv + Qvrg			
    
    Cmixven 	= (Qrm*Crm/PCmb + Qwm*Cwm/PCmb + Qvrg*Cvrg/PCvb + Qfat*Cfat/PCfb + Qliv*Cliv/PClib)/Qtot	
    Cart		= Clung/PClub																	
    
    dClung 	= (Qtot*(Cmixven-Clung/PClub) + Qalv*(Cair-Clung/PClub/PCba))/Vlung						
    dCrm 		= Qrm*(Cart-Crm/PCmb)/Vm															
    dCwm 		= Qwm*(Cart-Cwm/PCmb)/Vm															
    dCvrg 	= Qvrg*(Cart-Cvrg/PCvb)/Vvrg														
    dCfat 	= Qfat*(Cart-Cfat/PCfb)/Vfat														
    dCliv 	= (Qliv*(Cart-Cliv/PClib)-CLi*Cliv/PClib)/Vliv 										
    Cexh 		= QaQp* Cart/PCba+(1-QaQp)*Cair													
    
    Relupt 	= 100*Qalv*(Cair-Clung/PClub/PCba)/Qp/ppm*molvol										
    return(list(c(dClung, dCrm, dCwm, dCvrg, dCfat, dCliv), Cmixven = Cmixven, Cart = Cart, Cexh = Cexh, Relupt = Relupt))
  })
}

########################################################
# Initial concentrations:
yini <- c(	Clung = 0, 
           Crm = 0, 
           Cwm = 0,
           Cvrg = 0,
           Cfat = 0,
           Cliv = 0)

# Method to solve differential equations:
method 	= as.character(knime.in$"solvingmethod")
time_step = 0.1			

# Time course:
hours 	= 7	
start 	= 0			
stop 	= hours*60
times 	= seq(start, stop, time_step)

################ Input parameters  #############################

# Ventilation parameters (experimental):
Qp_work 		= knime.in$"Ventilation work"
Qp_rest 		= knime.in$"Ventilation rest"

# Exposure parameters:
molvol 		= knime.in$"Molar volume"
ppm 			= knime.in$"Exposure level"

# Physiological parameters:
BH 			= knime.in$"Body height" 
BW 			= knime.in$"Body weight"	

# Partition coefficents and clearance:
CLi 			= knime.in$"CLi"	
PCba 		= knime.in$"PCblood"	
PCfb 		= knime.in$"PCfat"	
PClub 		= knime.in$"PClung"	
PCmb 		= knime.in$"PCmuscle"	
PClib 		= knime.in$"PCliver"	
PCvb 		= knime.in$"PCvessel"

#################### Calculated Physiological parameters ##############

# Body composition
TBW 			= 12.86+0.1757*BH+0.331*BW		
FFM 			= TBW/0.72					
BV 			= FFM/1.1					

# Compartment volumes
Vfat 		= (BW-FFM)/0.92	
Vm 			= 0.344*BV		
Vlung 		= (0.00907+0.0193)*BV
Vliv 		= 0.0285*(TBW/0.72)	
Vbrain 		= 0.0256*BV		
Vkidney 		= 0.00532*BV			
Vother 		= 0.0103*BV		
Vvrg 		= Vbrain + Vkidney + Vother		

# Blood flows:
Qbrain 		= 0.57*Vbrain 					
Qother 		= 2.83*Vother					
Qrm 			= 0.0329*Vm						

QaQp 		= 0.89						

################### Input to list #############################

parameters <- c(
  Qp_work 	= Qp_work,
  Qp_rest 	= Qp_rest,
  molvol 		= molvol,
  ppm 		= ppm,
  BH 		= BH,
  BW 		= BW,
  CLi 		= CLi,
  PCba 		= PCba,
  PCfb 		= PCfb,	
  PClub		= PClub,
  PCmb		= PCmb,
  PClib		= PClib,
  PCvb 		= PCvb,
  TBW 		= TBW,
  FFM 		= FFM,
  BV 		= BV,
  Vfat 		= Vfat,
  Vm 		= Vm,
  Vlung		= Vlung,
  Vliv 		= Vliv,
  Vbrain		= Vbrain,
  Vkidney 	= Vkidney,
  Vother		= Vother,
  Vvrg 		= Vvrg, 
  Qbrain		= Qbrain,
  Qother		= Qother,
  Qrm 		= Qrm,
  QaQp		= QaQp)

############ Solving #######################################
solution <- ode(y = yini, times = times, func = pbpkEx, parms = parameters, method = method)
solution <- as.data.frame(solution)

############## Formatting ####################################
solution <- melt(solution, id.vars = "time", variable.name = "Compartment", value.name = "Concentration")

############## Going out #####################################
knime.out <- solution