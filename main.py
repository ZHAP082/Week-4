#This is an exact copy of the code from week 3. All you have to do is change the name of the file it reads at the begining to the galaxy whos graphs you wanna look at and change the rc value for that specific galaxy. I changed all the 1.87 for galaxy 1 to rc_value so you dont have to change each one you switch galaxies but you can use Ctrl F to find those.
rc_value = 1.87
import numpy as np
import matplotlib.pyplot as plt
############################################################################################################
#This bit of code is from week1 it makes an array for each label in Galaxy1.txt
############################################################################################################
f = open('Galaxy1 (4).txt','r')
list_make = f.readlines() 

All_data = []
Radias  = []
velocity = []
change_in_radias = []
change_in_velocity = []
Mass = []

for line in list_make:
    make_each_list = line.split('\t')
    All_data.append(make_each_list)

for i in All_data:
    Radias.append(i[0])
    velocity.append(i[1])
    change_in_radias.append(i[2])
    change_in_velocity.append(i[3])
    Mass.append(i[4])
    
Radias_Q = (Radias[1:])
velocity_Q = (velocity[1:])
change_in_radias_Q = (change_in_radias[1:])
change_in_velocity_Q = (change_in_velocity[1:])
Mass_data_n_Q = (Mass[1:])

Mass_data_Q = []

for i in Mass_data_n_Q:
    Just_data = i.replace('\n','')
    Mass_data_Q.append(Just_data)

Radias_data = []
velocity_data = []
change_in_radias_data = []
change_in_velocity_data = []
Mass_data = []

for i in Radias_Q:
    Radias_data.append(float(i))

for i in velocity_Q:
    velocity_data.append(float(i))

for i in change_in_radias_Q:
    change_in_radias_data.append(float(i))

for i in change_in_velocity_Q:
    change_in_velocity_data.append(float(i))

for i in Mass_data_Q:
    Mass_data.append(float(i))

Radias_data_array = np.array(Radias_data)
velocity_data_array = np.array(velocity_data)
change_in_radias_data_array = np.array(change_in_radias_data)
change_in_velocity_data_array = np.array(change_in_velocity_data)
Mass_data_array = np.array(Mass_data)
###############################################################################################################
###############################################################################################################


###############################################################################################################
#This is from week 2 it is plots for - v(r) for v observed, v visual, v sum. 
###############################################################################################################
plt.plot(Radias_data_array,velocity_data_array, label = 'Data')
plt.xlabel ('Radias / kpc')
plt.ylabel ('velocity / km/s')
plt.title ("Observed, calculated v's and uncertainty in optimum v against radius")

velocity_Visual_array = np.sqrt(((4.30*10**-6)*Mass_data_array)/Radias_data_array)

plt.plot(Radias_data_array,velocity_Visual_array, label = 'Visual mass')

Mass_DM_old = (4*np.pi*(100*10**6)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

Mass_sum = (Mass_data_array + Mass_DM_old)

velocity_sum_old = np.sqrt(((4.30*10**-6)*Mass_sum)/Radias_data_array) 

plt.plot(Radias_data_array,velocity_sum_old, label = 'Combined mass')
##################################################################################################################
##################################################################################################################

##################################################################################################################
#This is week 3. 
#1. Find X^2 for old p0 we were using - 100*10*6
#2. Iterate through range (guesses for optimum p0) and for each find X^2
#Then compare X^2 to old p0 X^2 and if X^2 is smaller means better that old p0 so that p0 and its X^2 added to better lists.
#Find smallest X^2 aka optimum p0.
#3. Using optimum p0 find optimum DM mass - find optimum sum mass - find optimum calculated v - use optimum v to plot optimum graph.
#4. This finds the uncertainty in optium density. It is the same code as 2. with a 2 at the end of each thing so no confusion and since we esentially do the same thing, but a diffrent if statment which finds all values of unceratniy. Then we find the mean of those uncertainties. This value is actually the effect of uncertainty on optimum density since we found densities for optimum X^2 + 1. Hence our actual uncertainty is the diffrence between that and our optimum density.
#5. Plot 2 graphs using optimum_density +- uncertainty to get +- velocity which use to graph. These are kinda like error bars. For my data turns put uncertainty is very small so in graph you have to zoom in alot to see uncertaity lines.
#6. Get fraction of visible mass in galaxy (visible / sum_mass) - many ways of doing this e.g get each fraction for each i and find mean of that. But i will get sum of all visible mass / sum of all combined_mass.
##################################################################################################################
#1.
Xsquared_values_old =  ((velocity_data_array - velocity_sum_old)**2)/(change_in_velocity_data_array **2)

Xsquared_old = np.sum(Xsquared_values_old)

better_guesses = []
better_Xsquared = []
##################################################################################################################
#2.

for density in np.arange((1*10**6), (200*10**6), 1000):

  Mass_density_DM_test = (4*np.pi*(density)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

  Cobmined_mass_test = (Mass_density_DM_test + Mass_data_array)

  v_model_test = np.sqrt(((4.30*10**-6)*Cobmined_mass_test)/Radias_data_array)

  Xsquared_values =  ((velocity_data_array - v_model_test)**2)/(change_in_velocity_data_array **2)

  X_squared = np.sum(Xsquared_values)

  if(X_squared < Xsquared_old):

    better_guesses.append(density)

    better_Xsquared.append(X_squared)
      
closest_0_X_squared = min(better_Xsquared, key=abs)

index_of_lowest_Xtwo = better_Xsquared.index(closest_0_X_squared)

Optimum_density =  (better_guesses[index_of_lowest_Xtwo])

print ('The optimum density is:')
print (Optimum_density)
#######################################################################################################################
#3.

Optimum_mass_DM = (4*np.pi*(Optimum_density)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

Optimum_Combined_mass = (Optimum_mass_DM + Mass_data_array)

v_optimum = np.sqrt(((4.30*10**-6)*Optimum_Combined_mass)/Radias_data_array)

plt.plot(Radias_data_array,v_optimum, label = 'Optimum combined mass')
########################################################################################################################
#4.MAYBE BREAKTHROUGH CHECK PARAGRAPH.

uncertainty_values = []

for density2 in np.arange((1*10**6), (200*10**6), 1000):

  Mass_density_DM_test2 = (4*np.pi*(density2)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

  Cobmined_mass_test2 = (Mass_density_DM_test2 + Mass_data_array)

  v_model_test2 = np.sqrt(((4.30*10**-6)*Cobmined_mass_test2)/Radias_data_array)

  Xsquared_values2 =  ((velocity_data_array - v_model_test2)**2)/(change_in_velocity_data_array **2)

  X_squared2 = np.sum(Xsquared_values2)

  if(int(X_squared2) == int((closest_0_X_squared + 1))):
    uncertainty_values.append(density2)
  
uncertainty_optium_density = (np.sum(uncertainty_values)) / (len(uncertainty_values))

actual_uncertainy_of_optimum_density = (uncertainty_optium_density - Optimum_density)

print('This is optimum_density due to uncertainty:') 
print(uncertainty_optium_density)
print('Hence the actual uncertainty is the diffrence between optimum_density due to uncertainty and optimum density: ')
print (actual_uncertainy_of_optimum_density)
########################################################################################################################
#5.

density_plus = (Optimum_density + actual_uncertainy_of_optimum_density)

density_minus = (Optimum_density - actual_uncertainy_of_optimum_density)

Mass_DM_plus = (4*np.pi*(density_plus)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

Mass_DM_minus = (4*np.pi*(density_minus)*(rc_value**2))*(Radias_data_array - (rc_value*np.arctan(Radias_data_array / rc_value)))

Mass_sum_plus = (Mass_DM_plus + Mass_data_array)

Mass_sum_minus = (Mass_DM_minus + Mass_data_array)

v_plus = np.sqrt(((4.30*10**-6)*Mass_sum_plus)/Radias_data_array)

v_minus = np.sqrt(((4.30*10**-6)*Mass_sum_minus)/Radias_data_array)

plt.plot(Radias_data_array,v_plus, label = '+ uncertainty')
plt.plot(Radias_data_array,v_minus, label = '- uncertainty')
##########################################################################################################################
#6.

sum_visible = np.sum(Mass_data_array)
sum_all = np.sum(Optimum_Combined_mass)

ratio_of_visible_mass_to_all_mass = ( sum_visible / sum_all )

print ('The fraction of visible mass to all mass in galaxy is :')
print (ratio_of_visible_mass_to_all_mass)

plt.legend(['Data', 'Visual mass', 'Combined mass', 'Optimum combined mass', '+ uncertainty', '- uncertainty'])

plt.show() # Have to put this here at the end since if you want to print someting it gas to be before pl.show.
############################################################################################################################