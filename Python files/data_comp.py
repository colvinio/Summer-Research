import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import datetime
import math
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



def are_holes(big_list, flag1, flag2, rad_cut):
     """
     Check for a gap of over some length in the center of the galaxy, plus at least 3 flags within
     """

     #Getting smallest y value
     metal_radii = big_list[1]
     smallest_y = min(metal_radii)
     
     #Counting flags within smallest_y
     count_flags = 0
     for radius in flag1[1]:
          if radius <= smallest_y:
               count_flags += 1
     for radius in flag2[1]:
          if radius <= smallest_y:
               count_flags += 1    

     #If there are no metal spaxids within and there are some number of flag spaxids within
          
     if .2 <= smallest_y < .3 and count_flags >= 6:
          return True
     elif .3 <= smallest_y < .45 and count_flags >= 8:
          return True
     elif .45 <= smallest_y and count_flags >= 10:
          return True
     else:
          return False



def get_line_vals(flag):
     """
     Returns [slope, intercept]
     """
     #create the line of best fit
     coefficients = np.polyfit(flag[1], flag[0], 1)
     slope = coefficients[0]
     intercept = coefficients[1]

     return [slope, intercept]



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



def galaxy_trash_can(flag_list1, flag_list2):
     """
     Checks if the galaxy is useable or not. Returns string "Trash" or "Not trash"
     """

     #make the new list so the originals lists are not edited
     double_checker = [[], []]
     for metal in flag_list1[0]:
          double_checker[0].append(metal)
     for metal in flag_list2[0]:
          double_checker[0].append(metal)
     for radius in flag_list1[1]:
          double_checker[1].append(radius)
     for radius in flag_list2[1]:
          double_checker[1].append(radius)

     # print()
     # print()
     # print(max(double_checker[0]) - min(double_checker[0]))
     # print(max(double_checker[1]) - min(double_checker[1]))
     
     #First: Are there enough spaxids? - takes 10 spaxids?
     if len(double_checker[0]) < 30:
          return "Too few spaxids"

     # #Second: Is the horizontal range large enough? Cannot be less than 1/2 of the vertical range?
     # #Horizontal radius can be no less than half of the vertical radius
     # horizontal_range = max(double_checker[1]) - min(double_checker[1])
     # vertical_range = max(double_checker[0]) - min(double_checker[0])
     # if vertical_range > 2 * horizontal_range:

     #New dynamic range way- min and max are factor of two apart at least
     if min(double_checker[1]) * 2 > max(double_checker[1]) and (max(double_checker[1]) - min(double_checker[1])) < .4:
          return "Bad range"
     
     else:
          return "Not trash"



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

      print(datetime.datetime.now().time())

      metal_types = ['Z94_smc', 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw', 'D16_smc', 'PP04_O3N2_smc', 'PP04_N2_smc', 'M13_N2_smc', 'M13_O3N2_smc', 'C17_O3N2_smc', 'C17_N2_smc', 'D16_mw', 'PP04_O3N2_mw', 'PP04_N2_mw', 'M13_N2_mw', 'M13_O3N2_mw', 'C17_O3N2_mw', 'C17_N2_mw']
      
      for type in metal_types:
          
          #Creating space for galaxies 'with holes' to have their IDs stored
          directory = f"ObjIDs, {type}"
          if not os.path.exists(f"{directory}"):
               os.makedirs(f"{directory}")

          #Make sure all files already in directory get deleted
          move_files_to_trash(directory)
          
          holy_objIDs = []
          unholy_IDs = []
          fewspaxs_IDS = []
          badrange_IDS = []

          #Opening SQL data
          data = open_file(f"/Users/colviniorio/Desktop/Research/Data/dr17 {type}.pkl")
          objIDs = data.keys()

          #limiter = 0
          
          #varius incrementers that will get used in print statements at the end
          num_thru = 0
          few_spaxs = 0
          bad_range = 0

          #For sorting galaxies later on
          no_holes = [[], []]
          yes_hole = [[], []]

          for ID in objIDs:
               
               vals_to_graph = organize_arrays(data[ID][0], data[ID][1], data[ID][2])

               S06_SF_vals = vals_to_graph.get('S06_SF')
               K03_SF_vals = vals_to_graph.get('K03_SF')
               K01_SF_vals = vals_to_graph.get('K01_SF')
               AGN_vals = vals_to_graph.get('AGN')

               #Safeguards any galaxy not having one of the flags
               if K01_SF_vals == None:
                    K01_SF_vals = [[], []]
               if AGN_vals == None:
                    AGN_vals = [[], []]
               if S06_SF_vals == None:
                    S06_SF_vals = [[], []]
               if K03_SF_vals == None:
                    K03_SF_vals = [[], []]
                    
               
               #Preparing the values to have avergaes found by getting rid of the nan values
               S06_SF, S06_SF_nans = destroy_nans(S06_SF_vals)
               K03_SF, K03_SF_nans = destroy_nans(K03_SF_vals)
               
               #Checking if the galaxy fits the criteria to be included in calculations --> if so, proceed with adding data
               if galaxy_trash_can(S06_SF, K03_SF) == "Not trash":
                    for metal in K03_SF[0]:
                         S06_SF[0].append(metal)
                    for radius in K03_SF[1]:
                         S06_SF[1].append(radius)

                    galaxy_data = get_line_vals(S06_SF)
                         
                    #Sort if galaxy has hole at center or not
                    hole_cut = .2
                    if are_holes(S06_SF, AGN_vals, K01_SF_vals, hole_cut) == True:
                         yes_hole[0].append(galaxy_data[0])
                         yes_hole[1].append(galaxy_data[1])
                         holy_objIDs.append(ID) #storing ID so we can graph and examine these graphs specifically
                    else:
                         no_holes[0].append(galaxy_data[0])
                         no_holes[1].append(galaxy_data[1])
                         unholy_IDs.append(ID) #storing ID so we can graph and examine these graphs specifically

               elif galaxy_trash_can(S06_SF, K03_SF) == "Too few spaxids":
                    fewspaxs_IDS.append(ID)

               else:
                    badrange_IDS.append(ID)

               num_thru += 1
               #print(f"{num_thru} out of {len(objIDs)} through.")
               #print()

               #  limiter += 1
               #  if limiter == 30:
               #       exit()

          #calculating averages to compare
          #try: #I only needed this when testing too few galaxies/too large a radius where there are no galaxies that 'have holes'
          average_yes = [sum(yes_hole[0]) / len(yes_hole[0]), sum(yes_hole[1]) / len(yes_hole[1])]
          average_no = [sum(no_holes[0]) / len(no_holes[0]), sum(no_holes[1]) / len(no_holes[1])]

          print()
          print(f"For {type}:")
          print(f"For {len(yes_hole[0])} galaxies with a hole in the center, average slope and intercept are {limit_significant_figures(average_yes[0], 4)} and {limit_significant_figures(average_yes[1], 4)}.")
          print(f"For the {len(no_holes[0])} other galaxies, average slope and intercept are {limit_significant_figures(average_no[0], 4)} and {limit_significant_figures(average_no[1], 4)}.")
          print(f"{len(fewspaxs_IDS)} galaxies were rejected for having 29 or fewer valid spaxids.")
          print(f"{len(badrange_IDS)} galaxies were rejected for their horizontal range being 'bad'.")
          print()
          
          #except:
          #    print('There were no galaxies that had holes of the specified length. I am too lazy to make this code still print stuff out right now.')
          
          #SAVE THE FILE
          pickle.dump(holy_objIDs, open(f"ObjIDs, {type}/ObjIDs of galaxies with holes.pkl", 'wb'))
          pickle.dump(unholy_IDs, open(f"ObjIDs, {type}/ObjIDs of galaxies without holes.pkl", 'wb'))
          pickle.dump(fewspaxs_IDS, open(f"ObjIDs, {type}/ObjIDs of galaxies, fewer than 30 spaxs.pkl", 'wb'))
          pickle.dump(badrange_IDS, open(f"ObjIDs, {type}/ObjIDs of galaxies, bad range.pkl", 'wb'))
          
      exit()


if __name__ == "__main__":
      main()