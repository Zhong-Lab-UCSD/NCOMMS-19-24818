#!/usr/bin/python3
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




import os
import sys
import glob
import math
from optparse import OptionParser







parameters = { 'InputFile':  None,
               'OutputPath': None,
               'chrSize':    None,
               'BaseName':   None,
               'reportUnit': 10000000,
               }







class ChrByChr:
    def __init__( self , parameters ):
        ##
        if os.path.exists( parameters['OutputPath'] ):
            pass
        else:
            os.makedirs( parameters['OutputPath'] )
        ##
        ##
        self.Read_chromosome_size( parameters )
        self.make_output_streams( parameters )
        self.Scan_table( parameters )
        self.close_output_streams( parameters )
        print(    '\n'    )


    def Read_chromosome_size( self , parameters ):
        self.chr____list = [  ]
        self.chr2size____dict = {  }
        ########
        print(    'Read chromosome size from: ' , parameters['chrSize']    )
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
            print(    '\t' , chromosome , '==>' , self.chr2size____dict[ chromosome ]    )
        ##
        self.chr____set = set( self.chr____list )
        print(    '\tChecking ==> ' , len( self.chr____set )    )
        

    def make_output_streams( self , parameters ):
        self.OutputFile____dict = {  }
        self.outputTEXT____dict = {  }
        ########
        print(    'Make output streams ......'    )
        for chrA in self.chr____list:
            for chrB in self.chr____list:
                chrA2chrBTuple = ( chrA , chrB )
                self.OutputFile____dict[ chrA2chrBTuple ]  =  parameters['OutputPath'] + parameters['BaseName'] +'.'+ chrA +'_'+ chrB +'.bedpe'
                print(    'Making stream ... ' , chrA2chrBTuple , '\t' , self.OutputFile____dict[ chrA2chrBTuple ]    )
                self.outputTEXT____dict[ chrA2chrBTuple ]  =  open( self.OutputFile____dict[ chrA2chrBTuple ] , 'w' )
            ##
        ##
        print( '' )


    def Scan_table( self , parameters ):
        count_line = 0            
        count_available_line = 0  
        count_available_line____dict = {  }  
        for chrA in self.chr____list:
            for chrB in self.chr____list:
                chrA2chrBTuple = ( chrA , chrB )
                count_available_line____dict[ chrA2chrBTuple ] = 0   
            ##
        ##
        ##
        print(    'Read from:  ' , parameters['InputFile']    )
        inputTEXT = open( parameters['InputFile'] , 'r' )    
        while True:
            try:
                line = next(inputTEXT)
                if line[0] == '#':
                    for chrA in self.chr____list:
                        for chrB in self.chr____list:
                            chrA2chrBTuple = ( chrA , chrB )
                            self.outputTEXT____dict[ chrA2chrBTuple ].write(  line  )
                        ##
                    ##
                elif line[0] != '#':
                    line2list = line.split('\n')[0].split('\r')[0].split('\t')
                    ##
                    count_line = count_line + 1    
                    ##
                    if count_line % parameters['reportUnit'] == 0:
                        print(    count_line , 'reads processed.'    )
                    else:
                        pass
                    ####
                    ####
                    chrA = line2list[0]
                    chrB = line2list[3]
                    ##
                    if ( chrA in self.chr____set ) and ( chrB in self.chr____set ):
                        count_available_line = count_available_line + 1  
                        chrA2chrBTuple = ( chrA , chrB )
                        self.outputTEXT____dict[ chrA2chrBTuple ].write(  line  )
                        count_available_line____dict[ chrA2chrBTuple ] = count_available_line____dict[ chrA2chrBTuple ] + 1    #### 计数
                    else:
                        pass
                    ##
                else:
                    pass
                ####
                if count_available_line == sum(count_available_line____dict.values()):
                    pass
                else:
                    print(    'Error ==> ' , count_available_line , sum(count_available_line____dict.values())    )
            except StopIteration:
                break
        inputTEXT.close()    
        print(    '\tChecking ==> Total lines: ' , count_line    )
        print(    '\tChecking ==> Total available lines: ' , count_available_line  ,   count_available_line == sum(count_available_line____dict.values())    )


    def close_output_streams( self , parameters ):
        print(    'Close output streams ......'    )
        for chrA in self.chr____list:
            for chrB in self.chr____list:
                chrA2chrBTuple = ( chrA , chrB )
                print(    'Closing stream ... ' , chrA2chrBTuple    )
                self.outputTEXT____dict[ chrA2chrBTuple ].close()
            ##
        ##
        print( '' )
        




if __name__ == '__main__':
    print(    '\n'    )


    
    usage='usage: python3 ChrByChr_split_of_iMARGI_reads.PY3.py [options]'
    parser = OptionParser( usage = 'python3  %prog  -i input_file    -o output_path    --name basename_of_outputfile    --chrSize chromosome_size_file' )
    parser.add_option('-i', '--input'     , dest='input_file'     , type='string', help='"/path/filename"  Input file.'            )
    parser.add_option('-o', '--output'    , dest='output_path'    , type='string', help='"/path/"  Output path.'                   )
    parser.add_option(      '--name'      , dest='base_name'      , type='string', help='Base name of outputfiles.'                )
    parser.add_option(      '--chrSize'   , dest='chromosome_size', type='string', help='"/path/filename"  Chromosome size file.'  )
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
    if options.base_name == None:
        parser.error('-h for help or provide the base name of output files!')
    else:
        pass
    ########
    parameters['InputFile']  = options.input_file
    parameters['OutputPath'] = options.output_path
    parameters['BaseName']   = options.base_name
    parameters['chrSize']    = options.chromosome_size
    
    

    ChrByChr( parameters )


    print(    '\nThe end.\n'    )













