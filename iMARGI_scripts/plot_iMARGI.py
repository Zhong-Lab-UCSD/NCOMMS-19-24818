#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser
import os
import numpy as np
import pandas as pd
import math

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

parameters = {'inputFile': None,
              'outputFile': None,
              'isGlobal': None,
              'chr_row': None,
              'chr_col': None,
              'bin_size': None,
              'select_part_map': None,
              'highlight_region': None,
              'data_type': None,
              'species': None,
              'custom_species': None,
              'my_colormap': None,
              'cutoff_type': None,
              'cutoff': None,
              'max_color': None,
              'my_dpi': None,
              'scale_factor': None,
              'plot_SE': None
              }

def save_matrix_space(input_matrix, output_filename):
    """
    Save a contact matrix in a txt file in a tab separated format. Columns are
    separated by tabs, rows are in different lines.
    Arguments:
        input_matrix (numpy matrix): input contact matrix to be saved
        output_filename (str): output file name in txt format
    Output:
        txt file containing the space separated data
    """
    import numpy as np
    rows = np.shape(input_matrix)[0]
    with open (output_filename, 'w') as f:
            for i in xrange(rows):
                row = [str(j) for j in input_matrix[i]]
                f.write(' '.join(row) + '\n')
                    

def convert_from_matrix_to_txt(input_file, output_txt):
    import numpy as np
    input_matrix = open(input_file, 'r')
    input_matrix_list = []
    while True:
        try:
            line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
            input_matrix_list.append(line2list[1:])
        except StopIteration:
            break
    input_matrix_array = np.array(input_matrix_list[1:]).astype(np.float)
    save_matrix_space(input_matrix_array,output_txt)
    
    
def convert_from_tab_to_txt(input_file, output_txt):
    import numpy as np
    input_matrix = open(input_file, 'r')
    input_matrix_list = []
    while True:
        try:
            line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
            input_matrix_list.append(line2list)
        except StopIteration:
            break
    input_matrix_array = np.array(input_matrix_list).astype(np.float)
    save_matrix_space(input_matrix_array,output_txt)

   
def extract_single_map(parameters):
    
    chr_row = parameters['chr_row']
    chr_col = parameters['chr_col']
    species = parameters['species']
    bin_size = int(parameters['bin_size'])
    data_type = parameters["data_type"]
    
    print("Extracting contact map chr" + chr_row + ' x chr' + chr_col + '...')
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
                chr_dim.append(int(line2list[1])/bin_size + 1) # + 1 is to take into account of the additional "final bin"
                d_chr_dim[line2list[0].replace('chr','')] = int(line2list[1])/bin_size + 1 # + 1 is to take into account of the additional "final bin"
            except StopIteration:
                break
            
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1

    input_file_list = map(str, parameters['inputFile'].strip('[]').split(','))
    index_file = 0
    for my_file in input_file_list:
        print("Loading " + my_file + '...')
        input_matrix = open(my_file, 'r')
        input_matrix_list = []
        index_file += 1
        while True:
            try:
                line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
                input_matrix_list.append(line2list[1:]) # to remove the row name (window)
            except StopIteration:
                break
        if index_file == 1:
            if parameters['scale_factor'] == None:
                input_matrix_array = np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
            else:
                input_matrix_array = np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
        else:
            if parameters['scale_factor'] == None:
                input_matrix_array = input_matrix_array + np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
            else:
                input_matrix_array = input_matrix_array + np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
        print("Done!")

    if chr_row == '1':
        row_start = 0
    else:
        row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_row)-1]]
    row_end = row_start + d_chr_dim[chr_row]
    
    if chr_col == '1':
        col_start = 0
    else:
        col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_col)-1]]
    col_end = col_start + d_chr_dim[chr_col]
    
    output_matrix = input_matrix_array[row_start:row_end,col_start:col_end]
    
    if chr_row == chr_col:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'iMARGI_' 'chr' + chr_row + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'iMARGI_' 'chr' + chr_row + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        save_matrix_space(output_matrix, my_filename)
    else:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        save_matrix_space(output_matrix, my_filename)
    print("Done!")
    return output_matrix


def plot_iMARGI_matrix (parameters):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    from matplotlib import rcParams
    
    #rcParams['font.family'] = 'arial'
    
    isGlobal = bool(int(parameters['isGlobal']))
    chr_row = parameters['chr_row']
    chr_col = parameters['chr_col']
    bin_size = int(parameters['bin_size'])
    data_type = parameters['data_type']
    species = parameters['species']
    my_dpi = int(parameters['my_dpi'])  
    
    if bin_size > 200000:
        grid_width = 2
    elif bin_size <= 200000 and bin_size > 100000:
        grid_width = 4
    elif bin_size <= 100000 and bin_size > 50000:
        grid_width = 8
    elif bin_size <= 50000:
        grid_width = 16
    
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
                chr_dim.append(int(line2list[1])/bin_size + 1) # + 1 is to take into account of the additional "final bin"
                d_chr_dim[line2list[0].replace('chr','')] = int(line2list[1])/bin_size + 1 # + 1 is to take into account of the additional "final bin"
            except StopIteration:
                break
            
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    d_chr_label_pos = {} # label position for chromosomes in the global matrix
    k = 0 # to consider the pixel occupied by the grid added after
    for i in chromosomes_list:
        d_chr_label_pos[i] = d_chr_dim_inc[i] - d_chr_dim[i]/2 + k
        k += grid_width # every grid line is grid_width
    
    label_pos = []
    label_name = []
    for i in chromosomes_list:
        label_pos.append(d_chr_label_pos[i])
        label_name.append('chr' + i)
    label_pos = np.array(label_pos)
    label_name = tuple(label_name)
    
    my_colormap = map(str, parameters['my_colormap'].strip('[]').split(','))
    if len(my_colormap) > 1:
        my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
    elif len(my_colormap) == 1:
        my_cmap = my_colormap[0]
    
    # Plot global heatmap
    if chr_row == None and chr_col == None:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'iMARGI_' + bin_size_str + '_' + data_type
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'iMARGI_' + bin_size_str + '_' + data_type
        
        input_file_list = map(str, parameters['inputFile'].strip('[]').split(','))
        index_file = 0
        for my_file in input_file_list:
            print("Loading " + my_file + '...')
            input_matrix = open(my_file, 'r')
            input_matrix_list = []
            index_file += 1
            while True:
                try:
                    line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
                    input_matrix_list.append(line2list[1:]) # to remove the row name (window)
                except StopIteration:
                    break
            if index_file == 1:
                if parameters['scale_factor'] == None:
                    matrix_data_full = np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
                else:
                    matrix_data_full = np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
            else:
                if parameters['scale_factor'] == None:
                    matrix_data_full = matrix_data_full + np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
                else:
                    matrix_data_full = matrix_data_full + np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
			
        print(str(np.sum(matrix_data_full)))
        print("Done!")
        print("Plotting " + parameters["inputFile"] + '...')
        
        # Adding grid to separate chromosomes
        k=0
        for i in chromosomes_list[:-1]:
            for j in range(grid_width):
                matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=1)
                matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=0)
            k += grid_width
        
        row = np.shape(matrix_data_full)[0]
        col = np.shape(matrix_data_full)[1]
        
        output_vect = np.reshape(matrix_data_full,row*col,1)
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
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == 'contact':
            perc = cutoff 
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == None:
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
            cbar = plt.colorbar()
            
        cbar.cmap.set_under('black')   
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.title(data_type + ' map (' + bin_size_str + ')', fontsize=12)
        plt.xticks(label_pos, label_name, rotation='vertical', fontsize = 6)
        plt.yticks(label_pos, label_name, fontsize = 6)
        plt.tick_params(axis='both', which='both', length=0)
        if parameters['outputFile'] == None:
            plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        else:
            plt.savefig(parameters['outputFile'], format = 'pdf', dpi=my_dpi)
        print("Done!")
    
    # Plot a single heatmap
    if (chr_row != None and chr_col == None) or (chr_row == None and chr_col != None):
        print("ERROR! Both the chromosomes have to be declared.")
        return
    if chr_row != None and chr_col != None:
        if isGlobal == True:
            matrix_data_full = extract_single_map(parameters) 
        else:
            input_file_list = map(str, parameters['inputFile'].strip('[]').split(','))
            index_file = 0
            for my_file in input_file_list:
                print("Loading " + my_file + '...')
                input_matrix = open(my_file, 'r')
                input_matrix_list = []
                index_file += 1
                while True:
                    try:
                        line2list = next(input_matrix).split('\n')[0].split('\r')[0].split('\t')
                        input_matrix_list.append(line2list[1:]) # to remove the row name (window)
                    except StopIteration:
                        break
                if index_file == 1:
                    if parameters['scale_factor'] == None:
                        matrix_data_full = np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
                    else:
                        matrix_data_full = np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
                else:
                    if parameters['scale_factor'] == None:
                        matrix_data_full = matrix_data_full + np.array(input_matrix_list[1:]).astype(np.float) # [1:] to skip the headline
                    else:
                        matrix_data_full = matrix_data_full + np.array(input_matrix_list[1:]).astype(np.float) * parameters["scale_factor"]
                print("Done!")

        chromosome_row = 'chr' + chr_row
        chromosome_col = 'chr' + chr_col
        
			
        print("Plotting " + chromosome_row + " x " + chromosome_col + " contact map...")
        
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'iMARGI_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'iMARGI_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type
        

        #### Highlighting region of the heatmap
        if parameters['highlight_region'] != None:
            print "Highlighting region..."
            region = map(int, parameters['highlight_region'].strip('[]').split(','))
            upper = int(math.floor(float(region[0])/bin_size)) - 1
            bottom = int(math.ceil(float(region[1])/bin_size))
            left = int(math.floor(float(region[2])/bin_size)) - 1
            right = int(math.ceil(float(region[3])/bin_size))
            # line at top of the region
            matrix_data_full[upper,left:right] = -1
            # line at bottom of the region
            matrix_data_full[bottom,left:right+1] = -1
            # line on left of the region
            matrix_data_full[upper:bottom,left] = -1
            # line on right of the region
            matrix_data_full[upper:bottom,right] = -1
        
        ### Loading SE track
        if parameters['plot_SE'] != None:
            temp = pd.read_csv(parameters['plot_SE'])
            se_track= list(temp["x"] - 1) # to 0-based coordinates
            for i in se_track:
                for j in xrange(np.shape(matrix_data_full)[1]):
                    if matrix_data_full[i,j] == 0:
                        matrix_data_full[i,j] = -1
                        
        #### Selecting a part of a single heatmap
        if parameters['select_part_map'] != None:
            select_part_map =  map(int, parameters['select_part_map'].strip('[]').split(','))
            if len(select_part_map) != 4:
                print("ERROR! Four coordinates should be inserted in this order: chr_row_start, chr_row_end, chr_col_start, chr_col_end.")
                return

            chr_row_coord = select_part_map[:2]
            chr_col_coord = select_part_map[2:]
            chr_row_bin = map(lambda x: float(x)/bin_size, chr_row_coord)
            chr_col_bin = map(lambda x: float(x)/bin_size, chr_col_coord)
            
            if chr_row_coord[0] >= chr_row_coord[1] or chr_col_coord[0] >= chr_col_coord[1]:
                print("ERROR! Start coordinate should be lower than end coordinate.")
                return
            if chr_row_bin[0] >= chr_row_bin[1] or chr_col_bin[0] >= chr_col_bin[1]:
                print("ERROR! Start coordinate should be much lower than the end coordinate given the bin size.")
                return
            if chr_row_bin[1] > d_chr_dim[chr_row]:
                print("ERROR! End coordinate of the chromosome on the rows should be lower than the chromosome size.")
                return
            if chr_col_bin[1] > d_chr_dim[chr_col]:
                print("ERROR! End coordinate of the chromosome on the cols should be lower than the chromosome size.")
                return
                
            matrix_data_full = matrix_data_full[int(math.floor(chr_row_bin[0])):int(math.ceil(chr_row_bin[1])),
                                                int(math.floor(chr_col_bin[0])):int(math.ceil(chr_col_bin[1]))]
        
        
        ########
        row = np.shape(matrix_data_full)[0]
        col = np.shape(matrix_data_full)[1]

        output_vect = np.reshape(matrix_data_full,row*col,1)
        non_zero = np.nonzero(output_vect)

        def format_e(n):
            a = '%e' % n
            return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]        
        
        plt.close("all")
        plt.gcf().subplots_adjust(left=0.15)
        plt.gcf().subplots_adjust(bottom=0.15)
        
        cutoff_type = parameters['cutoff_type']
        if parameters['cutoff_type'] != None:
            cutoff = float(parameters['cutoff'])
        max_color = parameters['max_color']
        
        if cutoff_type == 'perc':
            perc = np.percentile(output_vect[non_zero[0]],cutoff)
            if parameters['plot_SE'] == None:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
                cbar.cmap.set_under("#e6e6ff")
        elif cutoff_type == 'contact':
            perc = cutoff 
            if parameters['plot_SE'] == None:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
                cbar.cmap.set_under("#e6e6ff")
        elif cutoff_type == None:
            if parameters['plot_SE'] == None:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
                cbar = plt.colorbar()
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
                cbar = plt.colorbar()
                cbar.cmap.set_under("#e6e6ff")
        
        plt.title(data_type + ' map (' + bin_size_str + ')', fontsize=12)
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.ylabel(chromosome_row + ' coordinate (bp)', fontsize=10)
        plt.xlabel(chromosome_col + ' coordinate (bp)', fontsize=10)
        if parameters['select_part_map'] != None:
            ticks_row = (np.arange(0, row, row/4) * bin_size) + chr_row_coord[0]
            ticks_col = (np.arange(0, col, col/4) * bin_size) + chr_col_coord[0]
            format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
            format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
            plt.yticks(np.arange(0, row, row/4), format_ticks_row)
            plt.xticks(np.arange(0, col, col/4), format_ticks_col, rotation='vertical')
        else:
            ticks_row = (np.arange(0, row, row/4) * bin_size)
            ticks_col = (np.arange(0, col, col/4) * bin_size)
            format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
            format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
            plt.yticks(np.arange(0, row, row/4), format_ticks_row)
            plt.xticks(np.arange(0, col, col/4), format_ticks_col, rotation='vertical')
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
        if parameters['outputFile'] == None:
            plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        else:
            plt.savefig(parameters['outputFile'], format = 'pdf', dpi=my_dpi)
        
        print("Done!")
      

if __name__ == '__main__':
    usage='usage: python2.7 plot_iMARGI.py [options]'
    parser = OptionParser( usage = 'python2.7  %prog  -i inputFile -o outputFile')
    parser.add_option('-i', '--inputFile', dest='inputFile'  , type='string', help='"/path/filename.matrix"  Input file of iMARGI file matrix.')
    parser.add_option('-o', '--outputFile', dest='outputFile' , type='string', help='"/path/filename.pdf"  Output file of plot with pdf extension.')
    parser.add_option(      '--isGlobal', dest='isGlobal' , type='string', help='1 if inputFile is a global matrix (all-by-all chromosomes), 0 if it is a single contact matrix.')
    parser.add_option(      '--chr_row' , dest='chr_row' , type='str', help='Chromosome in the rows. Do not use it to plot the global matrix.')
    parser.add_option(      '--chr_col' , dest='chr_col' , type='str', help='Chromosome in the columns. Do not use it to plot the global matrix.')
    parser.add_option(      '--bin_size' , dest='bin_size' , type='str', help='Integer for the resolution of the contact matrices')
    parser.add_option(      '--select_part_map' , dest='select_part_map' , type='str', help='List of 4 elements as start_row_coordinate, end_row_coordinate, start_column_coordinate, end_column_coordinate (example: "[50000000,80000000,0,10000000]"')
    parser.add_option(      '--highlight_region' , dest='highlight_region' , type='str', help='List of of 4 elements as start_row_coordinate, end_row_coordinate, start_column_coordinate, end_column_coordinate (example: "[50000000,80000000,0,10000000]"')
    parser.add_option(      '--data_type' , dest='data_type' , type='str', help='Custom label to identify the type of data that you are plotting (example: observed, normalized, etc.).')
    parser.add_option(      '--species' , dest='species' , default='hg38', type='str', help='Species under analysis (hg38 or mm10). If a different species is used, insert the name here and provide "custom_species" parameter.')
    parser.add_option(      '--custom_species' , dest='custom_species' , type='str', help='"/path/filename.txt". Tab separated file for chromosomes of the custom species. Each row of the file is a chromosome: first column chromosome name with chr label (example:chrA); second column chromosome size.')
    parser.add_option(      '--my_colormap' , dest='my_colormap', default='[white,red]', type='str', help='List of colors to build the colormap from the lowest to the greatest number of contacts (at least two colors). A default colormap of matplotlib can be also chosen, in this case insert directly the colormap name. (Example: "[white,blue]" or "[#ffffff,#0000ff]" or "hot".')
    parser.add_option(      '--cutoff_type' , dest='cutoff_type', type='str', help='Upper cutoff type used for the number of contacts: it can be one of "percentile" (to use the percentile), "contact" (to use the number of contacts within a bin) or None to do not use any cutoff. Default: percentile.')
    parser.add_option(      '--cutoff' , dest='cutoff', default='95', type='str', help='If you are using a cutoff, insert here the cutoff (percentile as for example 99 for the 99th percentile or the number of contacts). Default: 95')
    parser.add_option(      '--max_color' , dest='max_color', default='#460000', type='str', help='If you use a cutoff, the color for the pixels above the cutoff.')
    parser.add_option(      '--my_dpi' , dest='my_dpi', default='2000', type='str', help='The resolution of the output figure in dpi (default 2000).')
    parser.add_option(      '--scale_factor' , dest='scale_factor', type='float', help='Scaling multiplicative factor to plot the heatmap.')
    parser.add_option(      '--plot_SE' , dest='plot_SE', type='string', help='File with the index of the contact matrix at the specific resolution to plot super enhancer stripes on the rows (to be calculated using R).')

    (options, args) = parser.parse_args( )
    ########
    if options.inputFile == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    ########
    parameters['inputFile'] = options.inputFile
    parameters['outputFile'] = options.outputFile
    parameters['isGlobal'] = options.isGlobal
    parameters['chr_row'] = options.chr_row
    parameters['chr_col'] = options.chr_col
    parameters['bin_size'] = options.bin_size
    parameters['select_part_map'] = options.select_part_map
    parameters['highlight_region'] = options.highlight_region
    parameters['data_type'] = options.data_type
    parameters['species'] = options.species
    parameters['custom_species'] = options.custom_species
    parameters['my_colormap'] = options.my_colormap
    parameters['cutoff_type'] = options.cutoff_type
    parameters['cutoff'] = options.cutoff
    parameters['max_color'] = options.max_color
    parameters['my_dpi'] = options.my_dpi
    parameters['scale_factor'] = options.scale_factor
    parameters['plot_SE'] = options.plot_SE

    plot_iMARGI_matrix(parameters)







