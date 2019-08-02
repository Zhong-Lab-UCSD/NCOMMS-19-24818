#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser
import os
import numpy as np

chromosomes = {'hg38':{'1':248956422,
                   '2':242193529,
                   '3':198295559,
                   '4':190214555,
                   '5':181538259,
                   '6':170805979,
                   '7':159345973,
                   '8':145138636,
                   '9':138394717,
                   '10':133797422,
                   '11':135086622,
                   '12':133275309,
                   '13':114364328,
                   '14':107043718,
                   '15':101991189,
                   '16':90338345,
                   '17':83257441,
                   '18':80373285,
                   '19':58617616,
                   '20':64444167,
                   '21':46709983,
                   '22':50818468,
                   'X':156040895,
                   'Y':57227415},
                   'mm10':{'1':195471971,
                   '2':182113224,
                   '3':160039680,
                   '4':156508116,
                   '5':151834684,
                   '6':149736546,
                   '7':145441459,
                   '8':129401213,
                   '9':124595110,
                   '10':130694993,
                   '11':122082543,
                   '12':120129022,
                   '13':120421639,
                   '14':124902244,
                   '15':104043685,
                   '16':98207768,
                   '17':94987271,
                   '18':90702639,
                   '19':61431566,
                   'X':171031299,
                   'Y':91744698}}  

parameters = {'inputFiles': None,
              'outputFile': None,
              'isGlobal': None,
              'chr_row': None,
              'chr_col': None,
              'time_points': None,
              'bin_size': None,
              #'select_part_map': None,
              'data_type': None,
              'species': None,
              'custom_species': None,
              'my_colormap': None,
              'cutoff_type': None,
              'cutoff': None,
              'max_color': None,
              'my_dpi': None,
              'topological_domains': None,
              'domain_color': None
              }


class plot_iMARGI_time:
    def __init__(self, parameters):
        self.plot_iMARGI_matrix(parameters)
        
    def plot_iMARGI_matrix (self, parameters):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np
        
        inputFiles = map(str, parameters['inputFiles'].strip('[]').split(','))
        isGlobal = bool(int(parameters['isGlobal']))
        chr_row = map(str, parameters['chr_row'].strip('[]').split(',')) # list of chromosomes in the rows
        chr_col = map(str, parameters['chr_col'].strip('[]').split(',')) # list of chromosomes in the columns
        time_points = map(str, parameters['time_points'].strip('[]').split(','))
        bin_size = int(parameters['bin_size'])
        data_type = parameters['data_type']
        species = parameters['species']
        my_dpi = int(parameters['my_dpi'])
        
        
        scaling_factors = {'Control': float(100000000)/float(65879147), 
                           'T3d_P6': float(100000000)/float(53291562), 
                           'T7d_P7': float(100000000)/float(50031654)}
        
        # Different width for the grid to separate contact maps at different resolutions in order to be well visualized in the plots
#        if bin_size > 200000:
#            grid_width = 2
#        elif bin_size <= 200000 and bin_size > 100000:
#            grid_width = 4
#        elif bin_size <= 100000 and bin_size > 50000:
#            grid_width = 8
#        elif bin_size <= 50000:
#            grid_width = 16
        grid_width = 2
        
        if species in chromosomes.keys():
            chromosomes_list = [str(i) for i in range(len(chromosomes[species]) - 1)[1:]] + ['X', 'Y']
            chr_dim = []
            for i in chromosomes_list:
                chr_dim.append(chromosomes[species][i]/bin_size + 1) # + 1 is to take into account of the additional "final bin"
            d_chr_dim = {}
            for i in chromosomes_list:
                d_chr_dim[i] = chromosomes[species][i]/bin_size + 1 # + 1 is to take into account of the additional "final bin"
        else:
            custom_chromosomes = open(parameters['custom_species'], 'r')
            chromosomes_list = []
            chr_dim = []
            d_chr_dim = {}
            while True:
                try:
                    line2list = next(custom_chromosomes).split('\n')[0].split('\t')
                    chromosomes_list.append(line2list[0].replace('chr',''))
                    chr_dim.append(int(line2list[1])/bin_size) # + 1 is to take into account of the additional "final bin"
                    d_chr_dim[line2list[0].replace('chr','')] = int(line2list[1])/bin_size # + 1 is to take into account of the additional "final bin"
                except StopIteration:
                    break
                
        d_chr_dim_inc = {}
        k=1
        for i in chromosomes_list:
            d_chr_dim_inc[i] = sum(chr_dim[:k])
            k+=1
        
        label_pos = []
        label_name = []
        last_pos = 0
        n = 0
        for chr_r,chr_c in zip(chr_row,chr_col):
            if n == 0:
                label_pos.append(d_chr_dim[chr_r]/2)
                last_pos += d_chr_dim[chr_r] + grid_width # to consider the grid
                n+=1
            else:
                label_pos.append(last_pos + d_chr_dim[chr_r]/2)
                last_pos += d_chr_dim[chr_r] + grid_width
                 
            if chr_r == chr_c:
                label_name.append('chr' + chr_r)
            else:
                label_name.append('chr' + chr_r + '-chr' + chr_c)
        label_pos = np.array(label_pos)
        label_name = tuple(label_name)
        
        my_colormap = map(str, parameters['my_colormap'].strip('[]').split(','))
        if len(my_colormap) > 1:
            my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
        elif len(my_colormap) == 1:
            my_cmap = my_colormap[0]
        
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'iMARGI_time' + bin_size_str + '_' + data_type
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'iMARGI_time' + bin_size_str + '_' + data_type
        
        # From global heatmap
        if isGlobal == True:
            # Load each global map at every time point into a dictionary
            inputFiles_dict = dict()
            for i in range(len(inputFiles)):
                print "Loading " + inputFiles[i] + '...'
                input_matrix = open(inputFiles[i], 'r')
                input_matrix_list = []
                while True:
                    try:
                        line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
                        input_matrix_list.append(line2list) # to remove the row name (window)
                    except StopIteration:
                        break
                matrix_data_full = np.array(input_matrix_list).astype(np.float)
                inputFiles_dict[time_points[i]] = matrix_data_full * scaling_factors[time_points[i]]
                print "Done!"
        
        time_steps = np.linspace(0,len(chr_row),11).astype(int).tolist() # to print percentage of completion to console


        print "(1/2) Building the matrix..."
        counter = 0 # to print percentage of completion to console
        init_counter = 0 # counter to initialize the output full matrix (if equal 1) or not
        n_col_list = [] # to save the uumber of columns in each row
        n_col_max = max([d_chr_dim[x] for x in chr_col]) * len(time_points) + 2*grid_width # to consider the grid
        for key, value in d_chr_dim.items():
            if value == max([d_chr_dim[x] for x in chr_col]):
                chr_col_max = key # bigger chromosomes in the columns
        output_full_matrix = np.zeros((1,n_col_max)) # initialize matrix so I can concatenate already a line where the chromosome in the columns is not the biggest
        for i,j in zip(chr_row, chr_col):
            init_counter += 1
            line_dict = dict() # to save the contact matrices per each line
            for k in time_points:
                # Extract the single map
                if isGlobal == True:
                    if i == '1':
                        row_start = 0
                    else:
                        row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(i)-1]]
                    row_end = row_start + d_chr_dim[i]
                    
                    if j == '1':
                        col_start = 0
                    else:
                        col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(j)-1]]
                    col_end = col_start + d_chr_dim[j]
                    
                    input_matrix_array = inputFiles_dict[k] # global map at the k time point
                    output_matrix = input_matrix_array[row_start:row_end,col_start:col_end]
                    #print str(output_matrix[25,51])
                
                else: # load each single file
                    input_matrix = open(inputFiles[time_points.index(k)] + 'chr' + i + '_chr' + j + '.normalized_byCol_byRow.matrix', 'r')
                    input_matrix_list = []
                    while True:
                        try:
                            line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
                            input_matrix_list.append(line2list[1:]) # to remove the row name (window)
                        except StopIteration:
                            break
                    output_matrix = np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
                    
                line_dict[k] = output_matrix # fill in the dictionary with each matrix per time point for this line
            
            # Build the line
            matrix_line = line_dict[time_points[0]] # initialize the line with the first matrix
            n_row = np.shape(matrix_line)[0]
            sep_col = np.zeros((n_row,grid_width))-1
            if j == chr_col_max: # I do not add the last sep_col since it's going till the full width of the heatmap
                for t in range(len(time_points))[1:]:
                    matrix_line = np.concatenate((matrix_line,sep_col,line_dict[time_points[t]]), axis=1)
            else: # for the other chromosomes that are shorter I add the last sep_col
                for t in range(len(time_points))[1:]:
                    if t != len(time_points)-1:
                        matrix_line = np.concatenate((matrix_line,sep_col,line_dict[time_points[t]]), axis=1)
                    else: # another sep_col is added at the end
                        matrix_line = np.concatenate((matrix_line,sep_col,line_dict[time_points[t]],sep_col), axis=1)
            
            # Attach the row to the output full matrix
            if init_counter == 1: # initialize the output full matrix
                n_col = np.shape(matrix_line)[1]
                diff_col = n_col_max - n_col
                matrix_line_full = np.concatenate((matrix_line,np.zeros((n_row,diff_col))), axis=1)
                output_full_matrix = np.concatenate((output_full_matrix,matrix_line_full), axis=0)
                n_col_list.append(n_col)
            else:
                n_col = np.shape(matrix_line)[1]
                diff_col = n_col_max - n_col
                matrix_line_full = np.concatenate((matrix_line,np.zeros((n_row,diff_col))), axis=1)
                if n_col < n_col_list[-1]: # the chromosome length is smaller than the previous chromosomes
                    sep_row = np.concatenate((np.zeros((grid_width,n_col_list[-1]))-1,np.zeros((grid_width,n_col_max-n_col_list[-1]))), axis=1)
                else:
                    sep_row = np.concatenate((np.zeros((grid_width,n_col))-1,np.zeros((grid_width,diff_col))), axis=1)
                n_col_list.append(n_col)
                output_full_matrix = np.concatenate((output_full_matrix,sep_row,matrix_line_full), axis=0)
            
            counter += 1
            for i in range(len(time_steps)):
                if counter == time_steps[i]:
                    print str(i*10) + '% completed.'
            

        print "(1/2) Done!"
    
        print "(2/2) Plotting..."
        
        output_full_matrix = output_full_matrix[1::]
        row = np.shape(output_full_matrix)[0]
        col = np.shape(output_full_matrix)[1]
        
        output_vect = np.reshape(output_full_matrix,row*col,1)
        negative_indexes = np.where(output_vect==-1)
        output_vect[negative_indexes] = 0
        non_zero = np.nonzero(output_vect)
        
        plt.close("all")
        
        cutoff_type = parameters['cutoff_type']
        if parameters['cutoff_type'] != None:
            cutoff = float(parameters['cutoff'])
        max_color = parameters['max_color']
        
        if cutoff_type == 'perc':
            perc = np.percentile(output_vect[non_zero[0]],cutoff)
            plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == 'contact':
            perc = cutoff 
            plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == None:
            plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmin=0)
            cbar = plt.colorbar()
        
        cbar.cmap.set_under('black')   
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.title(data_type + ' map SEs', fontsize=12)
        sample_pos = np.array([d_chr_dim[chr_col_max]/2, int(d_chr_dim[chr_col_max]*1.5), int(d_chr_dim[chr_col_max]*2.5)])
        plt.xticks(sample_pos, tuple(time_points), fontsize = 6)
        plt.yticks(label_pos, label_name, fontsize = 6)
        plt.tick_params(axis='both', which='both', length=0)
        if parameters['outputFile'] == None:
            plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        else:
            plt.savefig(parameters['outputFile'], format = 'pdf', dpi=my_dpi)

        print "(2/2) Done!"
        
#        # Plot a single heatmap
#        if (chr_row != None and chr_col == None) or (chr_row == None and chr_col != None):
#            print "ERROR! Both the chromosomes have to be declared."
#            return
#        if chr_row != None and chr_col != None:
#            if isGlobal == True:
#                matrix_data_full = self.extract_single_map(parameters)
#            else:
#                input_matrix = open(parameters['inputFile'], 'r')
#                input_matrix_list = []
#                while True:
#                    try:
#                        line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
#                        input_matrix_list.append(line2list[1:])
#                    except StopIteration:
#                        break
#                matrix_data_full = np.array(input_matrix_list[1:]).astype(np.float)
#
#            chromosome_row = 'chr' + chr_row
#            chromosome_col = 'chr' + chr_col
#            
#            print "Plotting " + chromosome_row + " x " + chromosome_col + " contact map..."
#            
#            if bin_size >= 1000000:
#                bin_size_str = str(bin_size/1000000) + 'mb'
#                my_filename = 'iMARGI_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type
#            elif bin_size < 1000000:
#                bin_size_str = str(bin_size/1000) + 'kb'
#                my_filename = 'iMARGI_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type
#            
#            #### Update matrix values to plot topological domains
#            if parameters['topological_domains'] != None:
#                if chr_row != chr_col:
#                    print "ERROR! To plot topological domains the matrix should be intrachromosomal."
#                    return
#                
#                domains_file = open(parameters['topological_domains'], 'r')
#                domains = []
#                while True:
#                    try:
#                        line2list = next(domains_file).split('\n')[0].split('\t')
#                        domains.append([int(x) for x in line2list])
#                    except StopIteration:
#                        break
#                my_filename = my_filename + '_domains'
#                diag_index = np.diag_indices(len(matrix_data_full))
#                for domain in domains:
#                    temp_start = domain[0]/bin_size
#                    temp_end = domain[1]/bin_size
#                    matrix_data_full[temp_start,temp_start:temp_end] = -1
#                    matrix_data_full[temp_start:temp_end,temp_end-1] = -1
#                    matrix_data_full[(diag_index[0][temp_start:temp_end],diag_index[1][temp_start:temp_end])] = -1
#            
#            #### Selecting a part of a single heatmap
#            if parameters['select_part_map'] != None:
#                select_part_map =  map(int, parameters['select_part_map'].strip('[]').split(','))
#                if len(select_part_map) != 4:
#                    print "ERROR! Four coordinates should be inserted in this order: chr_row_start, chr_row_end, chr_col_start, chr_col_end."
#                    return
#    
#                chr_row_coord = select_part_map[:2]
#                chr_col_coord = select_part_map[2:]
#                chr_row_bin = map(lambda x: x/bin_size, chr_row_coord)
#                chr_col_bin = map(lambda x: x/bin_size, chr_col_coord)
#                
#                if chr_row_coord[0] >= chr_row_coord[1] or chr_col_coord[0] >= chr_col_coord[1]:
#                    print "ERROR! Start coordinate should be lower than end coordinate."
#                    return
#                if chr_row_bin[0] >= chr_row_bin[1] or chr_col_bin[0] >= chr_col_bin[1]:
#                    print "ERROR! Start coordinate should be much lower than the end coordinate given the bin size."
#                    return
#                if chr_row_bin[1] > d_chr_dim[chr_row]:
#                    print "ERROR! End coordinate of the chromosome on the rows should be lower than the chromosome size."
#                    return
#                if chr_col_bin[1] > d_chr_dim[chr_col]:
#                    print "ERROR! End coordinate of the chromosome on the cols should be lower than the chromosome size."
#                    return
#                    
#                matrix_data_full = matrix_data_full[chr_row_bin[0]:chr_row_bin[1]+1,chr_col_bin[0]:chr_col_bin[1]+1]
#            
#            ########
#            row = np.shape(matrix_data_full)[0]
#            col = np.shape(matrix_data_full)[1]
#    
#            output_vect = np.reshape(matrix_data_full,row*col,1)
#            non_zero = np.nonzero(output_vect)
#    
#            def format_e(n):
#                a = '%e' % n
#                return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]        
#            
#            plt.close("all")
#            plt.gcf().subplots_adjust(left=0.15)
#            plt.gcf().subplots_adjust(bottom=0.15)
#            
#            cutoff_type = parameters['cutoff_type']
#            if parameters['cutoff_type'] != None:
#                cutoff = float(parameters['cutoff'])
#            max_color = parameters['max_color']
#            domain_color = parameters['domain_color']
#            
#            if cutoff_type == 'perc':
#                perc = np.percentile(output_vect[non_zero[0]],cutoff)
#                if parameters['topological_domains'] == None:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
#                    cbar = plt.colorbar(extend='max')
#                    cbar.cmap.set_over(max_color)
#                else:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
#                    cbar = plt.colorbar(extend='max')
#                    cbar.cmap.set_over(max_color)
#                    cbar.cmap.set_under(domain_color)
#            elif cutoff_type == 'contact':
#                perc = cutoff 
#                if parameters['topological_domains'] == None:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
#                    cbar = plt.colorbar(extend='max')
#                    cbar.cmap.set_over(max_color)
#                else:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
#                    cbar = plt.colorbar(extend='max')
#                    cbar.cmap.set_over(max_color)
#                    cbar.cmap.set_under(domain_color)
#            elif cutoff_type == None:
#                if parameters['topological_domains'] == None:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
#                    cbar = plt.colorbar()
#                else:
#                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
#                    cbar = plt.colorbar()
#                    cbar.cmap.set_under(domain_color)
#            
#            plt.title(data_type + ' map (' + bin_size_str + ')', fontsize=12)
#            cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
#            plt.ylabel(chromosome_row + ' coordinate (bp)', fontsize=10)
#            plt.xlabel(chromosome_col + ' coordinate (bp)', fontsize=10)
#            if parameters['select_part_map'] != None:
#                ticks_row = (np.arange(0, row, row/4) * bin_size) + chr_row_coord[0]
#                ticks_col = (np.arange(0, col, col/4) * bin_size) + chr_col_coord[0]
#                format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
#                format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
#                plt.yticks(np.arange(0, row, row/4), format_ticks_row)
#                plt.xticks(np.arange(0, col, col/4), format_ticks_col, rotation='vertical')
#            else:
#                ticks_row = (np.arange(0, row, row/4) * bin_size)
#                ticks_col = (np.arange(0, col, col/4) * bin_size)
#                format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
#                format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
#                plt.yticks(np.arange(0, row, row/4), format_ticks_row)
#                plt.xticks(np.arange(0, col, col/4), format_ticks_col, rotation='vertical')
#            plt.xticks(fontsize=8)
#            plt.yticks(fontsize=8)
#            plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
#            if parameters['outputFile'] == None:
#                plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
#            else:
#                plt.savefig(parameters['outputFile'], format = 'pdf', dpi=my_dpi)
#            
#            print "Done!"
      


if __name__ == '__main__':
    usage='usage: python2.7 plot_iMARGI_time.py [options]'
    parser = OptionParser( usage = 'python2.7  %prog  -i inputFiles -o outputFile --isGlobal --chr_row --chr_col --time_points --bin_size --data_type --species --custom_species --my_colormap --cutoff_type --cutoff --max_color --my_dpi --topological_domains --domain_color')
    parser.add_option('-i', '--inputFiles', dest='inputFiles'  , type='string', help='List with the input iMARGI global matrices, one per each condition to compare (each element as "/path/filename.matrix"). If isGlobal=False, insert the paths of the single contact matrices with the trailing slash "/" at the end.')
    parser.add_option('-o', '--outputFile', dest='outputFile' , type='string', help='"/path/filename.pdf"  Output file of plot with pdf extension.')
    parser.add_option(      '--isGlobal', dest='isGlobal' , type='string', help='1 if inputFiles contains global matrices, 0 if the inputFiles contains single contact matrices paths.')
    parser.add_option(      '--chr_row' , dest='chr_row' , type='str', help='List with the chromosomes in the rows (example: "[1,2,3]").')
    parser.add_option(      '--chr_col' , dest='chr_col' , type='str', help='List with the chromosomes in the columns (example: "[1,2,3]").')
    parser.add_option(      '--time_points' , dest='time_points' , type='str', help='List with the labels per each time point or condition (example: "[time1,time2]").')
    parser.add_option(      '--bin_size' , dest='bin_size' , type='str', help='Integer for the resolution of the contact matrices')
    parser.add_option(      '--data_type' , dest='data_type' , type='str', help='Custom label to identify the type of data that you are plotting (example: observed, normalized, etc.).')
    parser.add_option(      '--species' , dest='species' , default='hg38', type='str', help='Species under analysis (hg38 or mm10). If a different species is used, insert the name here and provide "custom_species" parameter.')
    parser.add_option(      '--custom_species' , dest='custom_species' , type='str', help='"/path/filename.txt". Tab separated file for chromosomes of the custom species. Each row of the file is a chromosome: first column chromosome name with chr label (example:chrA); second column chromosome size.')
    parser.add_option(      '--my_colormap' , dest='my_colormap', default='[white,red]', type='str', help='List of colors to build the colormap from the lowest to the greatest number of contacts (at least two colors). A default colormap of matplotlib can be also chosen, in this case insert directly the colormap name. (Example: "[white,blue]" or "[#ffffff,#0000ff]" or "hot".')
    parser.add_option(      '--cutoff_type' , dest='cutoff_type', default='perc', type='str', help='Upper cutoff type used for the number of contacts: it can be one of "percentile" (to use the percentile), "contact" (to use the number of contacts within a bin) or None to do not use any cutoff. Default: percentile.')
    parser.add_option(      '--cutoff' , dest='cutoff', default='95', type='str', help='If you are using a cutoff, insert here the cutoff (percentile as for example 99 for the 99th percentile or the number of contacts). Default: 95')
    parser.add_option(      '--max_color' , dest='max_color', default='#460000', type='str', help='If you use a cutoff, the color for the pixels above the cutoff.')
    parser.add_option(      '--my_dpi' , dest='my_dpi', default='2000', type='str', help='The resolution of the output figure in dpi (default 2000).')
    parser.add_option(      '--topological_domains' , dest='topological_domains', type='str', help='Not used for now.')
    parser.add_option(      '--domain_color' , dest='domain_color', default='#0000ff', type='str', help='Not used for now.')
    (options, args) = parser.parse_args( )
    ########
    if options.inputFiles == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    ########
    parameters['inputFiles'] = options.inputFiles
    parameters['outputFile'] = options.outputFile
    parameters['isGlobal'] = options.isGlobal
    parameters['chr_row'] = options.chr_row
    parameters['chr_col'] = options.chr_col
    parameters['time_points'] = options.time_points
    parameters['bin_size'] = options.bin_size
    parameters['data_type'] = options.data_type
    parameters['species'] = options.species
    parameters['custom_species'] = options.custom_species
    parameters['my_colormap'] = options.my_colormap
    parameters['cutoff_type'] = options.cutoff_type
    parameters['cutoff'] = options.cutoff
    parameters['max_color'] = options.max_color
    parameters['my_dpi'] = options.my_dpi
    parameters['topological_domains'] = options.topological_domains
    parameters['domain_color'] = options.domain_color

    plot_iMARGI_time(parameters)






