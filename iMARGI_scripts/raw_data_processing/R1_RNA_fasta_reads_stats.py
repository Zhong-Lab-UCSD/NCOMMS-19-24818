#!/usr/bin/python2.7
# -*- coding: utf-8 -*-



import sys
from optparse import OptionParser
import os

parameters = { 'InputFile':  None,
               'OutputFile': None,
               'Number':     None,
               }

def Fasta_stats( parameters ):
    print    'Filter fastq by marker ......'
    print    'Read from:  ' , parameters['InputFile']
    print    'Write down: ' , parameters['OutputFile']
    print    'Number:     ' , parameters['Number']
    ####
    code_dict = {  }
    count_line = 0               
    count_reads = 0              
    inputTEXT  = open( parameters['InputFile']  , 'r' )    
    while True:
        try:
            line = next(inputTEXT).split('\n')[0].split('\r')[0]
            count_line = count_line + 1    
            ####
            if count_line % 2 == 1:
                count_reads = count_reads + 1              
                ##
            elif count_line % 2 == 0:
                beginning = line[ 0 : parameters['Number'] ]
                ##
                if beginning not in code_dict.keys():
                    code_dict[  beginning  ] = 0
                else:
                    pass
                code_dict[  beginning  ] = code_dict[  beginning  ] + 1
                ##
        except StopIteration:
            break
    inputTEXT.close()    
    print    '\tChecking ==> Lines: ' , count_line 
    print    '\tChecking ==> Reads: ' , count_reads , '    ==> ' , count_reads , '* 2 =' , count_reads*2 , '    ==> ' , count_line == count_reads*2 , count_reads == sum( code_dict.values() ) 
    ########
    print    'Write down: ' , parameters['OutputFile']
    outputTEXT = open( parameters['OutputFile'] , 'w' )
    outputTEXT.write( 'Marker' +'\t'+ 'Count' +'\t'+ 'Percentage' +'\n' )
    for marker_count in sorted( code_dict.items() , key=lambda X:X[1] , reverse=True ):
        marker = marker_count[0]
        count  = str( marker_count[1] )
        percentage = str( float(count) / float(count_reads) )
        outputTEXT.write( marker +'\t'+ str(count) +'\t'+ str(percentage) +'\n' )
    outputTEXT.close()



if __name__ == '__main__':
    ####print    '\n'

    usage='usage: python R2_DNA_reads_filtering.py [options]'
    parser = OptionParser( usage = '%prog  -i /path/inputfile  -o /path/outputfile  -m nucleotide_marker  -n number_of_nucleotide' )
    parser.add_option('-i', '--input',  dest='input_file',                        type='string', help='"/path/filename"  Input file of fastq.')
    parser.add_option('-o', '--output', dest='output_file',                       type='string', help='"/path/filename"  Output file of statistics.')
    parser.add_option('-n', '--number', dest='number_of_nucleotide_as_marker',    type='int'   , help='Number of nucleotides at the beginning of the read as the marker for filter.')
    (options, args) = parser.parse_args( )
    ########
    if options.input_file == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    ########
    if options.output_file == None:
        parser.error('-h for help or provide the output file name!')
    else:
        pass
    ########
    if options.number_of_nucleotide_as_marker == None:
        parser.error('-h for help or provide the output file name!')
    else:
        pass
    ########
    parameters[ 'InputFile' ]  = options.input_file
    parameters[ 'OutputFile' ] = options.output_file
    parameters[ 'Number' ]     = options.number_of_nucleotide_as_marker
    

    Fasta_stats( parameters )

    
    ####print    '\nThe end.\n'










