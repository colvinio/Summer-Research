import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import math
import matplotlib.gridspec as gridspec
from send2trash import send2trash
import random
from scipy.interpolate import UnivariateSpline



def open_file(file_path):
     """
     Opens pickle file
     """
     
     # Open the pickle file in binary mode
     with open(file_path, "rb") as file:
          # Load the data from the pickle file
          data = pickle.load(file)
      
     return data 



def destroy_nans(big_list):
     """
     When fed a flag and corresponding lists, will get rid of all nans from main list
     The return will give back: big_list (without nans), nan_list (replacing nans with 0s so they can be graphed along axis without affecting the main slope and intercept)
     """

     nan_list = [[], []]
     to_delete = []

     #do this backwards so that in the next for loop, deleting first nans doesn't affect index of later nans
     for i in range(len(big_list[0]) - 1, -1, -1):
          #trying to get a boolean thats true for nans (which are floats?!?!)
          if big_list[0][i] != big_list[0][i]:
               to_delete.append(i)
      
     #to_delete is organized largest to smallest already
     for i in to_delete:
          nan_metal = big_list[0].pop(i)
          nan_radius = big_list[1].pop(i)
          
          #appends popped value to nan_list lists
          nan_list[0].append(nan_metal)
          nan_list[1].append(nan_radius)

     return big_list, nan_list



def sort_based_on_fir(line1, line2):
     """
     Sort and return two lists, based on the order that the first list should be sorted
     """

     sorted_list1, sorted_list2 = [], []

     extra_line1, extra_line2 = line1, line2

     for i in range(len(line1)):
        #Find index of minimum in line1
        min_index = extra_line1.index(min(extra_line1))
        
        #pop values from the duplicate lists at that index, store them in a new list
        smallest_x = extra_line1.pop(min_index)
        smallest_y = extra_line2.pop(min_index)
        
        sorted_list1.append(smallest_x)
        sorted_list2.append(smallest_y)

     return sorted_list1, sorted_list2



def find_cutoff_pt(derivatives, x_values):
     """
     Takes in derivatives and x_values, and returns singular x value where the derivative first becomes 3.5 times smaller than the initial value
     """

     first_slope = derivatives[0]

     for i in range(len(derivatives)):
          if first_slope < 0:
               if derivatives[i] > first_slope / 3.5:
                    #Find median slope pre-cut before returning values
                    med_slope = np.median(derivatives[:i])
                    med_xval = np.median(x_values[:i])
                    return x_values[i], med_slope, med_xval
          
     return "No cutoff pt", "No median index pre-cut", "No median x value pre-cut"



def percentile_cut_list(list, decimal):
     """
     Removes values outside of percentile specified (No return, edits lists themselves)
     """

     decimal_to_remove = (1 - decimal) / 2

     #Only keeping inner 95th percentile
     num_to_pop = round(decimal_to_remove * len(list))
     for i in range(num_to_pop):
        list.pop(0)
     for i in range(num_to_pop):
        list.pop(-1)

     

def sin_gal_graph(objID, S06_SF, K03_SF, K01_SF, AGN, directory, metal_type):
    """
    Graphs one galaxy
    Note: I HAVE TO FLIP THE 1 AND 0 INDICES BECAUSE I STORED METAL AND THEN RADIUS
    """

    #destroy nans
    S06_SF, S06_SF_nans = destroy_nans(S06_SF)
    K03_SF, K03_SF_nans = destroy_nans(K03_SF)

    if (len(S06_SF[0]) + len(K03_SF[0])) < 2:
         return

    else:
          #Height ratios between subplots (we want bottom subplot to be smaller)
          gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])  # Adjust height ratios as needed
          
          #first subplot
          ax1 = plt.subplot(gs[0])
          
          if len(S06_SF[0]) != 0:
               #GRAPHING S06_SF
               #plot non-nan values
               plt.scatter(S06_SF[1], S06_SF[0], c='red', label=f'S06_SF, {len(S06_SF_nans[0])} nulls', alpha = .3)

          if len(K03_SF[0]) != 0:
               #GRAPHING K03_SF
               #plot non-nan values
               plt.scatter(K03_SF[1], K03_SF[0], c='blue', marker="*", label=f'K03_SF, {len(K03_SF_nans[0])} nulls', alpha = .4)
          
          for metal in K03_SF[0]:
               S06_SF[0].append(metal)
          for radius in K03_SF[1]:
               S06_SF[1].append(radius)

        #   #Getting degree 2 poly of best fit
        #   coefficients = np.polyfit(S06_SF[1], S06_SF[0], 2)
        #   poly_func = np.poly1d(coefficients)
        #   x_values = np.linspace(min(S06_SF[1]), max(S06_SF[1]), 100)
        #   y_values = poly_func(x_values)

          sorted_xs, sorted_ys = sort_based_on_fir(S06_SF[1], S06_SF[0])
          
          #Only keeping the inner 95% of points
          percentile_cut_list(sorted_xs, .95), percentile_cut_list(sorted_ys, .95)
        
        # I don't think its working  #Doing stuff to try to smooth the spline out a little bit
        #   m = 2 * len(sorted_xs)
          
          #Getting spline of best fit (or trying to, at least)
          cs = UnivariateSpline(sorted_xs, sorted_ys) # Don't think this is doing anything: s=m+np.sqrt(2*m)
          x_values = np.linspace(min(sorted_xs), max(sorted_xs), 100)
          y_values = cs(x_values)
          
          plt.plot(x_values, y_values, color='grey', alpha=1, label=f"""Cubic spline of best fit""") #Label when using a poly of best fit: ax^2 + bx + c; a = {round(coefficients[0], 3)}; b = {round(coefficients[1], 3)}; c = {round(coefficients[2], 3)}

          #Finding instantaneous slopes at all x-values, then getting the median slope pre-cut
          array_derivatives = cs(sorted_xs, nu=1)
          cutoff, median_slope, median_x_val = find_cutoff_pt(array_derivatives, sorted_xs)

          #If theres a cutoff, plot vertical line, then plot line going through median pt with median slope
          if type(cutoff) is not str:
               #plot vertical line
               plt.plot([cutoff, cutoff], [min(sorted_ys), max(sorted_ys)], color='black', linestyle=':', alpha=.6, label=f'Cutoff @ x = {round(cutoff, 3)}')
               #Plot median line through median slope (There is so much typed in there because I need to find the exact y values for the two x points)
               plt.plot([min(sorted_xs), cutoff], [cs(median_x_val) - median_slope * (median_x_val - min(sorted_xs)), cs(median_x_val) + median_slope * (cutoff - median_x_val)], color='black', linestyle='--', linewidth=1.3, alpha=1, label=f'slope = {round(median_slope, 3)}')
          else:
               #Plot a line of best fit, but over all the data
               medianSlope = np.median(array_derivatives)
               medianPoint = np.median(x_values)
               medianYValue = cs(medianPoint)
               plt.plot([min(sorted_xs), max(sorted_xs)], [medianYValue + (min(sorted_xs) - medianPoint) * medianSlope, medianYValue + (max(sorted_xs) - medianPoint) * medianSlope], color='black', linestyle='--', linewidth=1.3, alpha=1, label=f'slope = {round(medianSlope, 3)}')

          #Wrapping up first plot display info
          plt.title(f"""ObjID {objID}, {metal_type}: Metallicity vs. Radius""")
          plt.ylabel("Z94_smc metallicity")
          plt.legend(fontsize='x-small')
          
          #Determing which of the two subplots has the larger x-axis value (for scaling)
          max_x = 0.1
          if len(K01_SF[1]) != 0:
               if max(K01_SF[1]) > max_x:
                    max_x = max(K01_SF[1])
          if len(AGN[1]) != 0:
               if max(AGN[1]) > max_x:
                    max_x = max(AGN[1])
          plt.xlim(0, max_x)

          # make tick labels invisible for first subplot
          plt.tick_params('x', labelbottom=False)
          
          #second subplot
          ax2 = plt.subplot(gs[1], sharex=ax1)
          plt.ylabel("Null flags")
          plt.xlabel("Radius from r_eff_nsa, units of sersic half light radii")

          #Plotting nulls from K01_SF and AGN
          #Get smallest value from S06_SF or K03_SF
          # y_val = min(S06_SF[0] + K03_SF[0])
          #Turn nans into y value
          for i in range(len(K01_SF[0])):
               K01_SF[0][i] = 0
          for i in range(len(AGN[0])):
               AGN[0][i] = 1
          plt.scatter(K01_SF[1], K01_SF[0], c='purple', alpha = .2)
          plt.scatter(AGN[1], AGN[0], c='green', marker="*", alpha = .35)

          #Various subplot editing things
          plt.yticks([0, 1], ['K01_SF', 'AGN'])
          plt.ylim(-1, 2)
          plt.xlim(0, max_x)
          
          # Adjust the spacing between subplots
          #plt.subplots_adjust(hspace=0.5)  # Increase/decrease as desired
          
          #Save and move on to next figure
          plt.savefig(f"{directory}/ObjID_{objID}, {metal_type}")
          plt.close()



def lines_vals_mxb(flag):
     """
     TAKES SWAPPED X AND Y: Function returns slope and intercept of line of best fit
     """

     #create the line of best fit
     coefficients = np.polyfit(flag[1], flag[0], 1)
     slope = coefficients[0]
     intercept = coefficients[1]

     return slope, intercept



def get_line(flag):
     """
     TAKES SWAPPED X AND Y: Function needs to have same horizontal (radius) bounds, create a line, and then give the two end points
     """

     #find horizontal bounds
     xmax = max(flag[1])
     xmin = min(flag[1])

     #create the line of best fit
     coefficients = np.polyfit(flag[1], flag[0], 1)
     slope = coefficients[0]
     intercept = coefficients[1]

     #find end points
     ymax = intercept + slope * xmax
     ymin = intercept + slope * xmin

     return [xmin, xmax], [ymin, ymax]



def organize_arrays(array1, array2, array3):
    """
    Function turns three jumbled arrays into a dictioary with 4 keys (flags) and the values being 2 lists of metallicities and radii
    """

    # Get unique flags from the third array
    flags = np.unique(array3)

    # Create an empty dictionary to store the organized data
    organized_data = {}

    # Iterate over the unique flags
    for flag in flags:
        # Get the indices where the flag matches in the third array
        indices = np.where(array3 == flag)[0]

        # Extract the corresponding values from array1 and array2
        values1 = array1[indices]
        values2 = array2[indices]

        # Append the flag and corresponding values to the organized data dictionary
        organized_data[flag] = [values1.tolist(), values2.tolist()]

    # Return the organized data
    return organized_data



def limit_significant_figures(number, significant_figures):
    if number == 0:
        return 0  # Special case for zero

    # Calculate the magnitude of the number
    magnitude = int(math.floor(math.log10(abs(number)))) 

    # Calculate the scaling factor based on the desired significant figures
    scaling_factor = 10 ** (significant_figures - magnitude - 1)

    # Scale the number, round it to the nearest integer, and scale it back
    rounded_number = round(number * scaling_factor) / scaling_factor

    return rounded_number



def move_files_to_trash(directory):
    """
    To be used to make sure that every time this program is run, only the most recent graphs remain
    """

    # Iterate over all the files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
       
        # Check if the current item is a file
        if os.path.isfile(file_path):
            
            # Move the file to the trash
            send2trash(file_path)
            print(f"Moved {file_path} to trash")



def main():

     #  print()
     #  print()
     #  print()
     #  print()
     #  print()

      metal_types = ['Z94_smc', 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw'] # , 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw'
      
      for metal_type in metal_types:
     
          #Opening the data
          data = open_file(f"/Users/colviniorio/Desktop/Research/Data/dr17 {metal_type}.pkl")
          
          # WHEN I WANTED TO GRAPH ALL 150 GALAXIES I INITIALLY PULLED
          #  #Make new directory
          #  directory = f"Graphs for first 150 galaxies"
          #  if not os.path.exists(f"{directory}"):
          #       os.makedirs(f"{directory}")
          #  objIDs = data.keys()
          
          direct_start = f"Graphs, {metal_type}"
          
          galID_types = ["galaxies with holes", "galaxies without holes"] #, "galaxies, fewer than 10 spaxs", "galaxies, bad range"
          
          for gallist in galID_types:

               # Accessing ID lists to pull correct info from data variable
               objIDs = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {metal_type}/ObjIDs of {gallist}.pkl")
               
               #Creating folder/accessing folder to put graphs into
               directory = f"{direct_start}/Graphs of {gallist}"
               if not os.path.exists(f"{directory}"):
                    os.makedirs(f"{directory}")

               #Make sure all files already in directory get deleted
               move_files_to_trash(directory)
               print()

               #Incrementor, used in print statements in the below for loop
               num_thru = 0
               
               random_vals = []

               num_of_graphs = 15
               
               for i in range(num_of_graphs): 
                    #Generate random number w/in length of objIDs list
                    random_vals.append(random.randint(0, len(objIDs) - 1))
                
               for i in random_vals:

                    vals_to_graph = organize_arrays(data[objIDs[i]][0], data[objIDs[i]][1], data[objIDs[i]][2])

                    S06_SF_vals = vals_to_graph.get('S06_SF')
                    K03_SF_vals = vals_to_graph.get('K03_SF')
                    K01_SF_vals = vals_to_graph.get('K01_SF')
                    AGN_vals = vals_to_graph.get('AGN')

                    #Safeguards any galaxy not having one of these two flags in the naming of the graph
                    if K01_SF_vals == None:
                         K01_SF_vals = [[], []]
                    if AGN_vals == None:
                         AGN_vals = [[], []]

                    #Should hopefully safeguard the graphing function if there are netiher of these flags?
                    if S06_SF_vals == None:
                         S06_SF_vals = [[], []]
                    if K03_SF_vals == None:
                         K03_SF_vals = [[], []]

                    #BE WARNED: below function EDITS THE LISTS THEMSELVES. the lists are no longer the same having been graphed
                    sin_gal_graph(objIDs[i], S06_SF_vals, K03_SF_vals, K01_SF_vals, AGN_vals, directory, metal_type)

                    num_thru += 1
                    print(f"For {metal_type}, {gallist}:")
                    print(f"{num_thru} out of {len(objIDs)} through (stopping at {num_of_graphs}).")
                    print()


      exit()


if __name__ == "__main__":
      main()