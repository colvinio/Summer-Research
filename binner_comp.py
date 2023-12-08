import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import math
import matplotlib.gridspec as gridspec
import random
import matplotlib.colors as colors
import matplotlib.cm as cm
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



def get_graphing_needs(data, objIDs, metalType):
    """
    Returns all the data that I need to graph all the plots
    return med_line_maker, density, xvalues, ylowvalues, yhighvalues, ylowestvalues, yhighestvalues, slopes, innerSplineYs, cutoffs
    """

    data_storer = []
    med_line_maker = [[], []]
    x_perc_line = []
    low_line = []
    high_line = []
    lowest_line = []
    highest_line = []
    density = []
    slopes = []
    cutoffs = []
    innerSplineYs = []
    
    for i in range(51):
        data_storer.append([[], []])

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

        # Combining data into one list
        for metal in K03_SF[0]:
            S06_SF[0].append(metal)
        for radius in K03_SF[1]:
            S06_SF[1].append(radius)
        
        #For line graphs- data storer is what splits data into bins
        for i in range(len(S06_SF[1])):
            if 0 <= S06_SF[1][i] <= 5:
                radius_mult = round(S06_SF[1][i] * 10)
                data_storer[radius_mult][0].append(S06_SF[1][i])
                data_storer[radius_mult][1].append(S06_SF[0][i])

        #CALCULATING SLOPE FOR HISTOGRAMS
        # Reduce spaxids to inner 95%
        x_values, y_values = sort_based_on_fir(S06_SF[1], S06_SF[0])
        percentile_cut_list(x_values, .95)
        percentile_cut_list(y_values, .95)
        
        # Get spline of best fit, then instantaneous slope at all points
        cs = UnivariateSpline(x_values, y_values)
        array_derivatives = cs(x_values, nu=1)

        #Also, add innermost y value of spline to the list
        innerSplineYs.append(cs(min(x_values)))

        # Check if there's a breakpoint
        cutoff, median_slope, median_x_val = find_cutoff_pt(array_derivatives, x_values)

        # Find and store median slope, either over inner 95% or from 2.5% to cutoff
        if type(cutoff) is not str:
            slopes.append(median_slope)
            cutoffs.append(cutoff)
        else:
            slopes.append(np.median(array_derivatives))

        num_thru += 1
        print(f"{num_thru} out of {len(objIDs)} through ({metalType}).")

    #RESUMING WITH LINE GRAPH - SOME CODE IS REDUNDANT, BUT IT'S A LITTLE TOO MUCH TO SORT THROUGH RIGHT NOW
    #Need to find x values for .025 and .975 percentiles
    x_value_list = []
    for i in range(len(data_storer)):
        for number in data_storer[i][0]:
            x_value_list.append(number)  
    x_value_list.sort()
    minimum_x_to_graph = np.percentile(x_value_list, 2.5)
    maximum_x_to_graph = np.percentile(x_value_list, 97.5)

    #Identify the bin the the min and max values fall within, then pop the bins that are outside of those bins
    reduced_data_storer = []
    index_low = math.floor(minimum_x_to_graph * 10)
    index_high = math.floor(maximum_x_to_graph * 10)
    for i in range(len(data_storer)):
        if i >= index_low and i <= index_high:
            reduced_data_storer.append([i, data_storer[i]])
    
    #Will have to plot the median line piece-wise to the the color to show up like i want it to
    for i in range(len(reduced_data_storer)):    
            #Find the median value, store it in final_line_maker
            median_metal = np.median(reduced_data_storer[i][1][1])
            med_line_maker[0].append(reduced_data_storer[i][0] / 10)
            med_line_maker[1].append(median_metal)

    #Find percentile of each bin
    for i in range(len(data_storer)):
        if len(data_storer[i][0]) != 0:
            #Storing x vals
            x_perc_line.append(i / 10)
            #storing inner 68% of vals
            low_line.append(np.percentile(data_storer[i][1], 16))
            high_line.append(np.percentile(data_storer[i][1], 84))
            #storing inner 95% of vals
            lowest_line.append(np.percentile(data_storer[i][1], 2.5))
            highest_line.append(np.percentile(data_storer[i][1], 97.5))
            #getting density
            density.append(len(data_storer[i][0]))
    
   #sorting data piece by piece so that I can change the colors
    xvalues = []
    ylowvalues = []
    yhighvalues = []
    ylowestvalues = []
    yhighestvalues = []
    
    for i in range(len(x_perc_line) - 1):
        xvalues.append([x_perc_line[i], x_perc_line[i + 1]])
        ylowvalues.append([low_line[i], low_line[i + 1]])
        yhighvalues.append([high_line[i], high_line[i + 1]])
        ylowestvalues.append([lowest_line[i], lowest_line[i + 1]])
        yhighestvalues.append([highest_line[i], highest_line[i + 1]])

    #HISTOGRAM STUFF AGAIN
    #Only keeping inner 95th percentile of slopes
    slopes.sort()
    percentile_cut_list(slopes, .95)

    #Cutting cutoff list to inner 95%
    cutoffs.sort()
    percentile_cut_list(cutoffs, .95)

    #Making density workable
    for i in range(len(density)):
        density[i] = math.log10(density[i])

    return med_line_maker, density, xvalues, ylowvalues, yhighvalues, ylowestvalues, yhighestvalues, slopes, innerSplineYs, cutoffs



def make_graph(namer, metal_type, medianLineMaker, spaxDensity, xValues, yLowVals, yHighVals, yLowestVals, yHighestVals):
    """
    This graphs
    """
    
    plt.title(f"{metal_type}, galaxies with {namer}", fontsize='small')
    plt.xlabel("Radius, from r_eff_nsa, units of sersic half light radii", fontsize='x-small')
    plt.ylabel(f"Metallicity from {metal_type}", fontsize='x-small')
    plt.plot(medianLineMaker[0], medianLineMaker[1], c='black', label = 'Median')

    #I continue trying valiantly to make this colorbar work
    #cmap = colors.Colormap('Reds', N=max(spaxDensity))
    cmap = cm.get_cmap('Reds') #, lut=max(spaxDensity)
    norm = colors.Normalize(vmin=min(spaxDensity), vmax=max(spaxDensity), clip=True)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, drawedges=False, location='bottom', pad=.2)
    cbar.ax.tick_params(labelsize='xx-small')
    cbar.set_label(label='Log Spaxids', fontsize='x-small')
    
   #Plot it piece by piece so that I can change the colors
    for i in range(len(xValues)):
        if i == 10:
            plt.fill_between(xValues[i], yLowVals[i], yHighVals[i], color=cmap(norm(spaxDensity[i])), alpha=.8, label = 'inner 68%') # cmap='rainbow'
            plt.fill_between(xValues[i], yLowestVals[i], yHighestVals[i], color=cmap(norm(spaxDensity[i])), alpha=.1, label = 'inner 95%') # cmap='rainbow'
        else:
            plt.fill_between(xValues[i], yLowVals[i], yHighVals[i], color=cmap(norm(spaxDensity[i])), alpha=.8)
            plt.fill_between(xValues[i], yLowestVals[i], yHighestVals[i], color=cmap(norm(spaxDensity[i])), alpha=.1)
    
    #Increasing number of ticks (autogenerated only shows at 0, 2 & 4)
    tick_positions = np.linspace(0, 5, 11)
    plt.xticks(tick_positions)

    plt.legend(fontsize='xx-small')
    plt.tight_layout()
    # Format the tick labels to display only whole numbers without decimals
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f'{int(value):.0f}' if value.is_integer() else f'{value:.1f}'))
    # Decrease the tick label size
    plt.tick_params(axis='both', which='both', labelsize='x-small')
    plt.grid()

    return



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



def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)



def make_histogram(data1, data2, title, xAxis, yAxes):
    """
    Holes data first, then no holes data. Making the histograms that don't have to do with the slope
    """

    #Make area under the curve 1
    #Find smallest and largest x val
    
    #First histogram
    plt.title(title, fontsize='xx-small')
    plt.xlabel(xAxis, fontsize='xx-small')
    plt.ylabel(yAxes, fontsize='xx-small')
    plt.tick_params(axis='both', which='both', labelsize='xx-small')
    plt.hist(data1, bins=30, density=True, edgecolor='blue', linewidth=1.5, label=f"Holes", histtype='step')

    #Overlayed histogram
    #plt.twinx()
    plt.hist(data2, bins=30, density=True, edgecolor="red", linewidth=1.5, alpha=.7, label=f"No holes", histtype='step')
    plt.tight_layout()
    plt.tick_params(axis='both', which='both', labelsize='xx-small')
    plt.legend(fontsize='xx-small', loc='upper left')

    return



def main():

     #  print()
     #  print()
     #  print()
     #  print()
     #  print()

     metal_types = ['Z94_smc', 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw'] # 
      
     for type in metal_types:
     
        #Opening the data
        data = open_file(f"/Users/colviniorio/Desktop/Research/Data/dr17 {type}.pkl")
    
        #Make new directory
        directory = f"Graphs, holes v no holes"
        if not os.path.exists(f"{directory}"):
            os.makedirs(f"{directory}")

        # pulling holes IDs
        objIDs_holes = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {type}/ObjIDs of galaxies with holes.pkl")
        holyMedLineMaker, holyDensity, holyXVals, holyYLow, holyYHigh, holyYLowest, holYHighest, holySlopes, holyInnerYs, holyCutoffs = get_graphing_needs(data, objIDs_holes, type)

        #pulling the non holes objIDs
        objIDs_noholes = open_file(f"/Users/colviniorio/Desktop/Research/ObjIDs, {type}/ObjIDs of galaxies without holes.pkl")
        unholyMedLineMaker, unolyDensity, unholyXVals, unholyYLow, unholyYHigh, unholyYLowest, unholYHighest, unholySlopes, unholyInnerYs, unholyCutoffs = get_graphing_needs(data, objIDs_noholes, type)
        
        #Initiate graph
        gs = gridspec.GridSpec(2, 6, height_ratios=[3, 2], hspace=.25, wspace=4)

        #Making colormap
        # custom_cmap = colors.LinearSegmentedColormap.from_list('custom', ['lightblue', 'blue', 'darkblue']
            
        #first subplot (holes graph)
        ax1 = plt.subplot(gs[0, :3])
        make_graph("holes", type, holyMedLineMaker, holyDensity, holyXVals, holyYLow, holyYHigh, holyYLowest, holYHighest)

        #second subplot (no holes graph)
        ax2 = plt.subplot(gs[0, 3:], sharey=ax1)
        make_graph('no holes', type, unholyMedLineMaker, unolyDensity, unholyXVals, unholyYLow, unholyYHigh, unholyYLowest, unholYHighest)
        
        #third subplot (slopes histogram)
        ax3 = plt.subplot(gs[1, :2])
        #To far in deep to make these more general... just going to make all new functions for the next two subplots
        make_histogram(holySlopes, unholySlopes, "Median slope of spline, inner 95th percentile", f"Metallicity, {type}", "Number of galaxies")

        #fourth subplot, innermost y value of spline
        ax4 = plt.subplot(gs[1, 2:4])
        make_histogram(holyInnerYs, unholyInnerYs, "Innermost y value of splines", f"Metallicity, {type}", "Number of galaxies")

        #Fifth subplot, cutoffs
        ax5 = plt.subplot(gs[1, 4:])
        make_histogram(holyCutoffs, unholyCutoffs, "Cutoff values, inner 95th percentile", "r_eff_rad", "Number of galaxies")
        
        #plt.subplots_adjust(hspace=.6)  # Increase/decrease as desired
        plt.savefig(f"{directory}/{type}- Gradient comp, holes vs. no hole.png")
        plt.close()

        print(f"Done with {type}!")
      
     exit()


if __name__ == "__main__":
      main()