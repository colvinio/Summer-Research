import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import math
import matplotlib.gridspec as gridspec
from send2trash import send2trash


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



def sin_gal_graph(objID, S06_SF, K03_SF, K01_SF, AGN, directory, type):
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
               plt.scatter(S06_SF[1], S06_SF[0], c='red', label=f'S06_SF, {len(S06_SF_nans[0])} nulls', alpha = .4)
               #line of best fit for non-nan data
               xData, yData = get_line(S06_SF)
               #plt.plot(xData, yData, linestyle='dotted', color = 'red', alpha = .8)

          if len(K03_SF[0]) != 0:
               #GRAPHING K03_SF
               #plot non-nan values
               plt.scatter(K03_SF[1], K03_SF[0], c='blue', marker="*", label=f'K03_SF, {len(K03_SF_nans[0])} nulls', alpha = .5)
               #line of best fit for non-nan data
               xData, yData = get_line(K03_SF)
               #plt.plot(xData, yData, linestyle='dashed', color = 'blue', alpha = .8)
          
          #Adding overall line of best fit
          for metal in K03_SF[0]:
               S06_SF[0].append(metal)
          for radius in K03_SF[1]:
               S06_SF[1].append(radius)
          xData, yData = get_line(S06_SF)
          slope, intercept = lines_vals_mxb(S06_SF)
          plt.plot(xData, yData, color = 'black', label=f'''slope={round(slope, 3)}, intercept={round(intercept, 3)}''')

          #Wrapping up first plot display info
          plt.title(f"""ObjID {objID}, {type}: Metallicity vs. Radius""")
          plt.ylabel("Z94_smc metallicity")
          plt.legend(fontsize='x-small')
          
          #Determing which of the two subplots has the larger x-axis value (for scaling)
          max_x = max(S06_SF[1])
          if len(K01_SF[1]) != 0:
               if max(K01_SF[1]) > max_x:
                    max_x = max(K01_SF[1])
          if len(AGN[1]) != 0:
               if max(AGN[1]) > max_x:
                    max_x = max(K01_SF[1])
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
          plt.scatter(K01_SF[1], K01_SF[0], c='purple', alpha = .4)
          plt.scatter(AGN[1], AGN[0], c='green', marker="*", alpha = .65)

          #Various subplot editing things
          plt.yticks([0, 1], ['K01_SF', 'AGN'])
          plt.ylim(-1, 2)
          plt.xlim(0, max_x)
          
          # Adjust the spacing between subplots
          #plt.subplots_adjust(hspace=0.5)  # Increase/decrease as desired
          
          #Save and move on to next figure
          plt.savefig(f"{directory}/ObjID_{objID}, {type}")
          plt.close()



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

      metal_types = ['Z94_smc', 'M91_smc'] # , 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw'
      
      for type in metal_types:
     
          #Opening the data
          data = open_file(f"/Users/colviniorio/Desktop/Research/Data/dr17 {type}.pkl")
          
          # WHEN I WANTED TO GRAPH ALL 150 GALAXIES I INITIALLY PULLED
          #  #Make new directory
          #  directory = f"Graphs for first 150 galaxies"
          #  if not os.path.exists(f"{directory}"):
          #       os.makedirs(f"{directory}")
          #  objIDs = data.keys()
          
          direct_start = f"Graphs, {type}"
          
          galID_types = ["galaxies with holes", "galaxies without holes", "galaxies, fewer than 10 spaxs", "galaxies, bad range"] #
          
          for gallist in galID_types:

               # Accessing ID lists to pull correct info from data variable
               objIDs = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {type}/ObjIDs of {gallist}.pkl")
               
               #Creating folder/accessing folder to put graphs into
               directory = f"{direct_start}/Graphs of {gallist}"
               if not os.path.exists(f"{directory}"):
                    os.makedirs(f"{directory}")

               #Make sure all files already in directory get deleted
               move_files_to_trash(directory)
               print()

               #Incrementor, used in print statements in the below for loop
               num_thru = 0
               
               for i in range(25):

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
                    #BEWARNED: below function EDITS THE LISTS THEMSELVES. the lists are no longer the same having been graphed
                    sin_gal_graph(objIDs[i], S06_SF_vals, K03_SF_vals, K01_SF_vals, AGN_vals, directory, type)

                    num_thru += 1
                    print(f"For {type}, {gallist}:")
                    print(f"{num_thru} out of {len(objIDs)} through (stopping at 25).")
                    print()


      exit()


if __name__ == "__main__":
      main()