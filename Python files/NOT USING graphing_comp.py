import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import math
import matplotlib.gridspec as gridspec
import random


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



def lines_vals_mxb(flag):
     """
     Function returns slope and intercept of line of best fit
     """

     #create the line of best fit
     coefficients = np.polyfit(flag[1], flag[0], 1)
     slope = coefficients[0]
     intercept = coefficients[1]

     return slope, intercept



def get_line(flag):
     """
     Function needs to have same horizontal (radius) bounds, create a line, and then give the two end points
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



def make_graph(data, objIDs, namer, metal_type):
    """
    Using this to make the graphs with the bazillion lines
    """

    #list used to store all values, [[slopes], [intercepts], [xmins], [xmaxes], [ymins], [ymaxes]]
    data_storer = [[], [], [], [], [], []]
    
    plt.title(f"{metal_type} gradients, {len(objIDs)} galaxies with {namer}", fontsize='small')
    plt.xlabel("Radius, from r_eff_nsa, units of sersic half light radii", fontsize='x-small')
    plt.ylabel(f"Metallicity from {metal_type}", fontsize='x-small')

    num_thru = 0
      
    for ID in objIDs:

        vals_to_graph = organize_arrays(data[ID][0], data[ID][1], data[ID][2])

        S06_SF_vals = vals_to_graph.get('S06_SF')
        K03_SF_vals = vals_to_graph.get('K03_SF')

        #Should hopefully safeguard the graphing function if there are netiher of these flags?
        if S06_SF_vals == None:
            S06_SF_vals = [[], []]
        if K03_SF_vals == None:
            K03_SF_vals = [[], []]
        
        S06_SF, S06_SF_nans = destroy_nans(S06_SF_vals)
        K03_SF, K03_SF_nans = destroy_nans(K03_SF_vals)

        for metal in K03_SF[0]:
            S06_SF[0].append(metal)
        for radius in K03_SF[1]:
            S06_SF[1].append(radius)

        #calculate slope, intercept, and x-bounds of galaxy
        xvals, yvals = get_line(S06_SF)
        slope, intercept = lines_vals_mxb(S06_SF)

        #store these values
        data_storer[0].append(slope)
        data_storer[1].append(intercept)
        data_storer[2].append(xvals[0])
        data_storer[3].append(xvals[1])
        data_storer[4].append(yvals[0])
        data_storer[5].append(yvals[1])

        #changing colors slightly
        colTuple = (.5 + (random.randrange(0, 50, 1) - 25)/100, .5 + (random.randrange(0, 50, 1) - 25)/100, .5 + (random.randrange(0, 50, 1) - 25)/100)

        #Plot the line of best fit for this specific galaxy on the graph
        plt.plot(xvals, yvals, color = colTuple, alpha = .3)

        num_thru += 1
        print(f"{num_thru} out of {len(objIDs)} through ({metal_type}).")


    #Plot average line of best fit on graph
    mean_slope = np.mean(data_storer[0])
    mean_intercept = np.mean(data_storer[1])

    #generate line values
    x_line = np.linspace(min(data_storer[2]), max(data_storer[3]), 100)
    y_line = mean_slope * x_line + mean_intercept

    #Check y limits (only need to check min, because should have a negative slope)
    if min(y_line) < min(data_storer[4]):
        x_max = (min(data_storer[4]) - mean_intercept) / mean_slope
        x_line = np.linspace(min(data_storer[2]), x_max, 100)
        y_line = mean_slope * x_line + mean_intercept

    #plot line
    plt.plot(x_line, y_line, color = 'black', linestyle = 'dotted', label = f'''mean slope and intercept
    slope={limit_significant_figures(mean_slope, 4)} , intercept={limit_significant_figures(mean_intercept, 4)}''')

    #plot median line of best fit on graph
    median_slope = np.median(data_storer[0])
    median_intercept = np.median(data_storer[1])

    #generate line values
    x_line = np.linspace(min(data_storer[2]), max(data_storer[3]), 100)
    y_line = median_slope * x_line + median_intercept

    #Check y limits (only need to check min, because should have a negative slope)
    if min(y_line) < min(data_storer[4]):
        x_max = (min(data_storer[4]) - median_intercept) / median_slope
        x_line = np.linspace(min(data_storer[2]), x_max, 100)
        y_line = median_slope * x_line + median_intercept

    #plot line
    plt.plot(x_line, y_line, color = 'black', linestyle = 'dashed', label = f'''median slope and intercept
    slope={limit_significant_figures(median_slope, 4)}, intercept={limit_significant_figures(median_intercept, 4)}''')

    #save plot
    plt.legend(fontsize='xx-small')
    plt.tight_layout()

    return



def make_histogram(data, objIDs, metal_type):
     """
     Makes the histogram
     """

     print("Making histogram...")
     
     plt.title(f"{metal_type} gradient, 95th percentile", fontsize='small')
     plt.xlabel(f"{metal_type} gradient slope", fontsize='x-small')
     plt.ylabel("Number of galaxies", fontsize='x-small')

     num_thru = 0
     slopes = []
      
     for ID in objIDs:

        vals_to_graph = organize_arrays(data[ID][0], data[ID][1], data[ID][2])

        S06_SF_vals = vals_to_graph.get('S06_SF')
        K03_SF_vals = vals_to_graph.get('K03_SF')

        #Should hopefully safeguard the graphing function if there are netiher of these flags?
        if S06_SF_vals == None:
            S06_SF_vals = [[], []]
        if K03_SF_vals == None:
            K03_SF_vals = [[], []]
        
        S06_SF, S06_SF_nans = destroy_nans(S06_SF_vals)
        K03_SF, K03_SF_nans = destroy_nans(K03_SF_vals)

        for metal in K03_SF[0]:
            S06_SF[0].append(metal)
        for radius in K03_SF[1]:
            S06_SF[1].append(radius)

        #calculate slope
        slope, intercept = lines_vals_mxb(S06_SF)

        #store slope
        slopes.append(slope)

     #Organizing slopes from largest to smallest
     slopes.sort()

     #Only keeping inner 95th percentile
     num_to_pop = round(.025 * len(slopes))
     for i in range(num_to_pop):
        slopes.pop(0)
     for i in range(num_to_pop):
        slopes.pop(-1)
     
     plt.hist(slopes, bins=30, edgecolor='black')
     plt.tight_layout()



def main():

     #  print()
     #  print()
     #  print()
     #  print()
     #  print()

     metal_types = ['Z94_smc', 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw'] # , 'M91_mw', 'KK04_mw', 'KE08_mw'
      
     for type in metal_types:
     
        #Opening the data
        data = open_file(f"/Users/colviniorio/Desktop/Research/Data/dr17 {type}.pkl")
        
    
        #Make new directory
        directory = f"Graphs, holes v no holes"
        if not os.path.exists(f"{directory}"):
            os.makedirs(f"{directory}")

        # puuling holes IDs
        objIDs_holes = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {type}/ObjIDs of galaxies with holes.pkl")

        #pulling the non holes objIDs
        objIDs_noholes = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {type}/ObjIDs of galaxies without holes.pkl")
        
        #Initiate graph
        gs = gridspec.GridSpec(2, 2)
            
        #first subplot (holes graph)
        ax1 = plt.subplot(gs[0])
        make_graph(data, objIDs_holes, "holes", type)
        
        #second subplot (no holes graph)
        ax2 = plt.subplot(gs[1])
        make_graph(data, objIDs_noholes, 'no holes', type)

        #third subplot (holes histogram)
        ax3 = plt.subplot(gs[2])
        make_histogram(data, objIDs_holes, type)

        #fourth subplot (no holes histogram)
        ax4 = plt.subplot(gs[3])
        make_histogram(data, objIDs_noholes, type)
        
        plt.subplots_adjust(hspace=.6)  # Increase/decrease as desired
        plt.savefig(f"{directory}/{type}- Gradient comp, holes vs. no hole.png")
        plt.close()

        print(f"Done with {type}!")
      
     exit()


if __name__ == "__main__":
      main()