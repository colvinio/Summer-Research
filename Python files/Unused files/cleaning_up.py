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


def get_graphing_needs(data, objIDs, metalType):
    """
    Returns all the data that I need to graph all the plots
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
    for i in range(len(x_perc_line) - 1):
        xvalues = [x_perc_line[i], x_perc_line[i + 1],]
        ylowvalues = [low_line[i], low_line[i + 1]]
        yhighvalues = [high_line[i], high_line[i + 1]]
        ylowestvalues = [lowest_line[i], lowest_line[i + 1]]
        yhighestvalues = [highest_line[i], highest_line[i + 1]]

    #HISTOGRAM STUFF AGAIN
    #Organizing slopes from largest to smallest
    slopes.sort()

    #Only keeping inner 95th percentile
    num_to_pop = round(.025 * len(slopes))
    for i in range(num_to_pop):
        slopes.pop(0)
    for i in range(num_to_pop):
        slopes.pop(-1)

    return med_line_maker, density, xvalues, ylowvalues, yhighvalues, ylowestvalues, yhighestvalues, slopes, innerSplineYs, cutoffs



def make_graph(namer, metal_type, cmap, medianLineMaker, spaxDensity, xValues, yLowVals, yHighVals, yLowestVals, yHighestVals):
    """
    This graphs
    """
    
    plt.title(f"{metal_type}, {len(objIDs)} galaxies with {namer}", fontsize='small')
    plt.xlabel("Radius, from r_eff_nsa, units of sersic half light radii", fontsize='x-small')
    plt.ylabel(f"Metallicity from {metal_type}", fontsize='x-small')

    plt.plot(medianLineMaker[0], medianLineMaker[1], c='black', label = 'Median')
    
   #Plot it piece by piece so that I can change the colors
    for i in range(xValues):
        if i == 10:
            plt.fill_between(xValues, yLowVals, yHighVals, color=cmap(spaxDensity[i]), alpha=.8, label = 'inner 68%') # cmap='rainbow'
            plt.fill_between(xValues, yLowestVals, yHighestVals, color=cmap(spaxDensity[i]), alpha=.1, label = 'inner 95%') # cmap='rainbow'
        else:
            plt.fill_between(xValues, yLowVals, yHighVals, color=cmap(spaxDensity[i]), alpha=.8)
            plt.fill_between(xValues, yLowestVals, yHighestVals, color=cmap(spaxDensity[i]), alpha=.1)
    
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




#SAVE THESE FUNCTIONS INCASE I MESS SOMETHING UP DRASTICALLY

def make_slope_histogram(data, objIDs, metal_type):
     """
     Makes the histogram
     """

     print("Making histogram...")

     num_thru = 0
     slopes = []
     cutoffs = []
     innerSplineYs = []
      
     number_through = 0
     
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

        #CALCULATING SLOPE
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

        number_through += 1
        print(f'{number_through} galaxies out of {len(objIDs)} galaxies through the {metal_type} histogram')

     #Organizing slopes from largest to smallest
     slopes.sort()

     #Only keeping inner 95th percentile
     num_to_pop = round(.025 * len(slopes))
     for i in range(num_to_pop):
        slopes.pop(0)
     for i in range(num_to_pop):
        slopes.pop(-1)
     
     plt.title(f"{metal_type} gradient, inner 95th percentile", fontsize='xx-small')
     plt.xlabel(f"{metal_type} gradient slope", fontsize='xx-small')
     plt.ylabel("Number of galaxies", fontsize='xx-small')
     plt.hist(slopes, bins=30, edgecolor='black', color="blue", label=f"Holes", histtype='step', fill=True)
     plt.legend(fontsize='xx-small', loc='upper left')
     # Decrease the tick label size
     plt.tick_params(axis='both', which='both', labelsize='xx-small')

     return innerSplineYs, cutoffs



def overlay_slope_histogram(data, objIDs, metal_type):
     """
     Overlays the other data on top
     """

     plt.twinx()
     
     num_thru = 0
     slopes = []
     cutoffs = []
     innerSplineYs = []
      
     number_through = 0
     
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

        #CALCULATING SLOPE
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

        number_through += 1
        print(f'{number_through} galaxies out of {len(objIDs)} galaxies through the {metal_type} histogram')

     #Organizing slopes from largest to smallest
     slopes.sort()

     #Only keeping inner 95th percentile
     num_to_pop = round(.025 * len(slopes))
     for i in range(num_to_pop):
        slopes.pop(0)
     for i in range(num_to_pop):
        slopes.pop(-1)
     
     plt.hist(slopes, bins=30, edgecolor='black', color="red", alpha=.5, label=f"No holes", histtype='step', fill=True)
     plt.tight_layout()
     plt.tick_params(axis='both', which='both', labelsize='xx-small')
     plt.legend(fontsize='xx-small')

     return innerSplineYs, cutoffs



def make_graph(data, objIDs, namer, metal_type, cmap):
    """
    This graph should make 50 bins and then find the median values at 0, .1, .2, etc up to 5. Then connect the points. For error, this uses percentile bounds.
    """

    data_storer = []
    med_line_maker = [[], []]
    x_perc_line = []
    low_line = []
    high_line = []
    lowest_line = []
    highest_line = []
    density = []
    
    for i in range(51):
        data_storer.append([[], []])
    
    plt.title(f"{metal_type}, {len(objIDs)} galaxies with {namer}", fontsize='small')
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

        # Combining data into one list
        for metal in K03_SF[0]:
            S06_SF[0].append(metal)
        for radius in K03_SF[1]:
            S06_SF[1].append(radius)
        
        for i in range(len(S06_SF[1])):
            if 0 <= S06_SF[1][i] <= 5:
                radius_mult = round(S06_SF[1][i] * 10)
                data_storer[radius_mult][0].append(S06_SF[1][i])
                data_storer[radius_mult][1].append(S06_SF[0][i])

        num_thru += 1
        print(f"{num_thru} out of {len(objIDs)} through ({metal_type}).")

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

    #Plot the median line
    plt.plot(med_line_maker[0], med_line_maker[1], c='black', label = 'Median')

    # #Finding density
    # for i in range(len(data_storer)):
    #     if len(data_storer[i][0]) < 10:
    #         density.append(0)
    #     else:
    #         density.append(round(math.log(len(data_storer[i][0]))))
    
    for i in range(len(data_storer)):
        density.append(len(data_storer[i][0]))

    #log_density = np.log10(density)

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
    
   #Plot it piece by piece so that I can change the colors
    for i in range(len(x_perc_line) - 1):
        
        xvalues = [x_perc_line[i], x_perc_line[i + 1],]
        ylowvalues = [low_line[i], low_line[i + 1]]
        yhighvalues = [high_line[i], high_line[i + 1]]
        ylowestvalues = [lowest_line[i], lowest_line[i + 1]]
        yhighestvalues = [highest_line[i], highest_line[i + 1]]

        if i == 10:
            plt.fill_between(xvalues, ylowvalues, yhighvalues, color=cmap(density[i]), alpha=.8, label = 'inner 68%') # cmap='rainbow'
            plt.fill_between(xvalues, ylowestvalues, yhighestvalues, color=cmap(density[i]), alpha=.1, label = 'inner 95%') # cmap='rainbow'
        else:
            plt.fill_between(xvalues, ylowvalues, yhighvalues, color=cmap(density[i]), alpha=.8)
            plt.fill_between(xvalues, ylowestvalues, yhighestvalues, color=cmap(density[i]), alpha=.1)
    
    #Increasing number of ticks (autogenerated only shows at 0, 2 & 4)
    tick_positions = np.linspace(0, 5, 11)
    plt.xticks(tick_positions)
    
    #save plot
    plt.legend(fontsize='xx-small')
    plt.tight_layout()
    # Format the tick labels to display only whole numbers without decimals
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f'{int(value):.0f}' if value.is_integer() else f'{value:.1f}'))
    # Decrease the tick label size
    plt.tick_params(axis='both', which='both', labelsize='x-small')
    plt.grid()

    return density