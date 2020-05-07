#!/usr/bin/python2.7
# -*- coding: utf-8 -*-





'''
hg38
'''
'''
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
chr11	135086622
chr12	133275309
chr13	114364328
chr14	107043718
chr15	101991189
chr16	90338345
chr17	83257441
chr18	80373285
chr19	58617616
chr20	64444167
chr21	46709983
chr22	50818468
chrX	156040895
chrY	57227415
'''




import sys
from optparse import OptionParser
import os




parameters = { 'InputFile':  None,
               'OutputPath': None,
               'chrSize':    None,
               'cutoff':     None,
               'resolution': None,
               'chromA':     None,
               'chromB':     None,
               'reportUnit': 10000000,
               }








class Make_matrix:
    def __init__(  self  ,  parameters  ):
        ##
        if os.path.exists( parameters['OutputPath'] ):
            pass
        else:
            os.makedirs( parameters['OutputPath'] )
        ##
        self.Read_chromosome_size( parameters )
        self.Make_window( parameters )
        self.Read_MARGI( parameters )
        self.Write_down( parameters )
        

    def Read_chromosome_size( self , parameters ):
        self.chr____list = [  ]
        self.chr2size____dict = {  }
        ########
        print    'Read chromosome size from: ' , parameters['chrSize']
        inputTEXT = open( parameters['chrSize'] , 'r' )    
        while True:
            try:
                line2list = next(inputTEXT).split('\n')[0].split('\r')[0].split('\t')
                chromosome = line2list[0]
                chrSize    = int( line2list[1] )
                ##
                self.chr____list.append( chromosome )
                self.chr2size____dict[ chromosome ] = chrSize
            except StopIteration:
                break
        inputTEXT.close()    
        for chromosome in self.chr____list:
            print    '\t' , chromosome , '==>' , self.chr2size____dict[ chromosome ]
        ##
        self.chr____set = set( self.chr____list )
        print    '\tChecking ==> ' , len( self.chr____set )
        

    def Make_window( self , parameters ):
        self.chrAwindow_list = [  ]
        self.chrBwindow_list = [  ]
        self.chrAwindow_chrBwindow_count____2D_dict = {  }
        ########
        print    'Make windows for: ' , parameters['chromA']
        chromSize = self.chr2size____dict[ parameters['chromA'] ]
        for start in range(  1 , chromSize + 1 , parameters['resolution']  ):
            end = min( start + parameters['resolution'] - 1    ,    chromSize )
            self.chrAwindow_list.append(    ( parameters['chromA'] , start , end )    )
        print    '\tChecking ==> self.chrAwindow_list: ' , len( self.chrAwindow_list )
        ##
        ##
        print    'Make windows for: ' , parameters['chromB']
        chromSize = self.chr2size____dict[ parameters['chromB'] ]
        for start in range(  1 , chromSize + 1 , parameters['resolution']  ):
            end = min( start + parameters['resolution'] - 1    ,    chromSize )
            self.chrBwindow_list.append(    ( parameters['chromB'] , start , end )    )
        print    '\tChecking ==> self.chrBwindow_list: ' , len( self.chrBwindow_list )
        ##
        ##
        print    'Build window to window 2D dict ...... '
        count_line = 0    
        for chrAwindow in self.chrAwindow_list:
            self.chrAwindow_chrBwindow_count____2D_dict[ chrAwindow ] = {  }
            count_line = count_line + 1    
            sys.stdout.write('\r        Line:'+str(count_line) )
            for chrBwindow in self.chrBwindow_list:
                self.chrAwindow_chrBwindow_count____2D_dict[ chrAwindow ][ chrBwindow ] = 0
            sys.stdout.flush()
        print    '\tChecking ==> Lines: ' , count_line
        print    '\tChecking ==> self.chrAwindow_chrBwindow_count____2D_dict: ' , len( self.chrAwindow_chrBwindow_count____2D_dict.keys() )


    def Read_MARGI( self , parameters ):
        print    'Read from: ' , parameters['InputFile']
        count_line = 0              
        count_available_line = 0      
        count_unavailable_line = 0    
        ########
        inputTEXT = open( parameters['InputFile'] , 'r' )    
        while True:
            try:
                line = next(inputTEXT).split('\n')[0].split('\r')[0]
                if line[0] != '#':
                    line2list = line.split('\t')
                    count_line = count_line + 1    
                    ##
                    if count_line % parameters['reportUnit'] == 0:
                        print    count_line , 'reads processed.'
                    else:
                        pass
                    ##
                    RNA_chromosome = line2list[0]
                    RNA_start  = int(line2list[1]) + 1    #### coordinate in 1-based system
                    RNA_end    = int(line2list[2])        #### coordinate in 1-based system
                    ##
                    DNA_chromosome = line2list[3]
                    DNA_start  = int(line2list[4]) + 1    #### coordinate in 1-based system
                    DNA_end    = int(line2list[5])        #### coordinate in 1-based system
                    ####
                    if line2list[-1] == 'Inter-chromosome':
                        d = 10000000000000000
                    else:
                        d = int( line2list[-1] )
                    ####
                    if ( RNA_chromosome == parameters['chromA'] ) and ( DNA_chromosome == parameters['chromB'] ) and ( d > parameters['cutoff'] ):
                        RNA_chromSize = self.chr2size____dict[ RNA_chromosome ]
                        DNA_chromSize = self.chr2size____dict[ DNA_chromosome ]
                        ########
                        count_available_line = count_available_line + 1      
                        ##
                        RNA_midPoint = ( RNA_start + RNA_end ) // 2
                        RNA_window_start = RNA_midPoint // parameters['resolution'] * parameters['resolution'] + 1            #### coordinate in 1-based system
                        RNA_window_end   = RNA_window_start + parameters['resolution'] - 1                                    #### coordinate in 1-based system
                        RNA_window = ( RNA_chromosome  ,  max( RNA_window_start , 1 )  ,  min( RNA_window_end , RNA_chromSize ) )
                        ##
                        DNA_midPoint = ( DNA_start + DNA_end ) // 2
                        DNA_window_start = DNA_midPoint // parameters['resolution'] * parameters['resolution'] + 1            #### coordinate in 1-based system
                        DNA_window_end   = DNA_window_start + parameters['resolution'] - 1                                    #### coordinate in 1-based system
                        DNA_window = ( DNA_chromosome  ,  max( DNA_window_start , 1 )  ,  min( DNA_window_end , DNA_chromSize ) )
                        ##
                        self.chrAwindow_chrBwindow_count____2D_dict[  RNA_window  ][  DNA_window  ] = self.chrAwindow_chrBwindow_count____2D_dict[  RNA_window  ][  DNA_window  ] + 1    ########    recording    ########
                        ##
                    else:
                        count_unavailable_line = count_unavailable_line + 1    
                    ####
                    ####
                    if count_available_line + count_unavailable_line == count_line:
                        pass
                    else:
                        print    'Error A'
                    ##
                    ##
                else:
                    pass
            except StopIteration:
                break
        inputTEXT.close()    
        print    '\tChecking ==> Lines: ' , count_line
        print    '\tChecking ==> Available lines: ' , count_available_line
        print    '\tChecking ==> Total: ' , count_available_line + count_unavailable_line , count_available_line + count_unavailable_line == count_line


    def Write_down( self , parameters ):
        OutputFile_matrix = parameters['OutputPath'] + parameters['chromA'] +'_'+ parameters['chromB'] + '.matrix'
        OutputFile_cdt    = parameters['OutputPath'] + parameters['chromA'] +'_'+ parameters['chromB'] + '.cdt'
        print    'Write down matrix:  ' , OutputFile_matrix
        print    'Write down heatmap: ' , OutputFile_cdt
        count_line = 0    
        ##
        outputTEXT_matrix = open( OutputFile_matrix , 'w' )
        outputTEXT_matrix.write(    'Matrix' +'\t'+ ('\t').join( [ window[0]+':'+str(window[1])+'-'+str(window[2]) for window in self.chrBwindow_list ] ) + '\n'    )    #### chromosome+':'+str(start)+'-'+str(end)
        ##
        outputTEXT_cdt = open( OutputFile_cdt , 'w' )
        outputTEXT_cdt.write(    ('\t').join( ['Matrix' ,'NAME','GWEIGHT']  +  [ window[0]+':'+str(window[1])+'-'+str(window[2]) for window in self.chrBwindow_list ] ) +'\n'    )    #### chromosome+':'+str(start)+'-'+str(end)
        outputTEXT_cdt.write(    ('\t').join( ['EWEIGHT',  ''  ,    ''   ]  +  ['1']*len(self.chrBwindow_list)                                                        ) +'\n'    )    #### chromosome+':'+str(start)+'-'+str(end)
        ##
        for RNA_window in self.chrAwindow_list:    #### chromosome+':'+str(start)+'-'+str(end)
            count_line = count_line + 1    
            sys.stdout.write('\r        Line:'+str(count_line) )
            ##
            row_head = RNA_window[0]+':'+str(RNA_window[1])+'-'+str(RNA_window[2])
            written_data_list = [  str( self.chrAwindow_chrBwindow_count____2D_dict[ RNA_window ][ DNA_window ] )  for  DNA_window  in  self.chrBwindow_list  ]
            ##
            outputTEXT_matrix.write(    row_head +'\t'+ ('\t').join( written_data_list ) + '\n'    )
            outputTEXT_cdt.write(    ('\t').join( [ row_head , row_head , '1' ]  +  written_data_list    ) +'\n'    )
            ##
            del row_head
            del written_data_list
            ##
            sys.stdout.flush()
        ##
        outputTEXT_matrix.close()
        outputTEXT_cdt.close()
        print    '\tChecking ==> Lines: ' , count_line

            
            

if __name__ == '__main__':
    print    '\n'



    usage='usage: python2.7 IntraChrom_RNAvsDNA_matrix.PY27.py [options]'
    parser = OptionParser( usage = 'python2.7  %prog  -i input_file    -o output_path    --chrSize chromosome_size_file    -c cutoff    -r resolution' )
    parser.add_option('-i', '--input'     , dest='input_file'     , type='string', help='"/path/filename"  Input file of iMARGI bedpe file.'  )
    parser.add_option('-o', '--output'    , dest='output_path'    , type='string', help='"/path/"  Output path of matrix.'                    )
    parser.add_option(      '--chrSize'   , dest='chromosome_size', type='string', help='"/path/filename"  Chromosome size file.'             )
    parser.add_option('-c', '--cutoff'    , dest='distance_cutoff', type='int'   , help='Threshold to distinguish between distal and proximal interactions for intra-chromosome read pairs.' )
    parser.add_option('-r', '--resolution', dest='resolution'     , type='int'   , help='Resolution of output matrix.' )
    parser.add_option(      '--rna_chr'   , dest='rna_chr'        , type='string', help='chromosome of RNA end.'       )
    parser.add_option(      '--dna_chr'   , dest='dna_chr'        , type='string', help='chromosome of DNA end.'       )
    (options, args) = parser.parse_args( )
    ########
    if options.input_file == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    ########
    if options.output_path == None:
        parser.error('-h for help or provide the output path!')
    else:
        pass
    ########
    if options.chromosome_size == None:
        parser.error('-h for help or provide the chromosome size file name!')
    else:
        pass
    ########
    if options.distance_cutoff == None:
        parser.error('-h for help or provide the cutoff threshold!')
    else:
        pass
    ########
    if options.resolution == None:
        parser.error('-h for help or provide the resolution of matrix!')
    else:
        pass
    ########
    if options.rna_chr == None:
        parser.error('-h for help or provide the RNA chromosome!')
    else:
        pass
    ########
    if options.dna_chr == None:
        parser.error('-h for help or provide the DNA chromosome!')
    else:
        pass
    ########
    parameters['InputFile']  = options.input_file
    parameters['OutputPath'] = options.output_path
    parameters['chrSize']    = options.chromosome_size
    parameters['cutoff']     = options.distance_cutoff
    parameters['resolution'] = options.resolution
    parameters['chromA']     = options.rna_chr
    parameters['chromB']     = options.dna_chr
    

    Make_matrix( parameters )


    print    '\nThe end.\n'













