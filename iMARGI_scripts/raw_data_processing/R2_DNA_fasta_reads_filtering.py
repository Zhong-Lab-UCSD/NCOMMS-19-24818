#!/usr/bin/python2.7
# -*- coding: utf-8 -*-



import sys
from optparse import OptionParser
import os

parameters = { 'InputFile':  None,
               'OutputFile': None,
               'Marker':     None,
               'Number':     None,
               }

def Fasta_filter( parameters ):
    print    'Filter fastq by marker ......'
    print    'Read from:  ' , parameters['InputFile']
    print    'Write down: ' , parameters['OutputFile']
    print    'Marker:     ' , parameters['Marker']
    print    'Number:     ' , parameters['Number']
    ####
    Output_discard = parameters['OutputFile'] + '.discarded'
    ####
    code_dict = {  }
    count_line = 0               
    count_reads = 0              
    count_available_line = 0     
    count_available_reads = 0    
    count_unavailable_line = 0   
    count_unavailable_reads = 0  
    inputTEXT  = open( parameters['InputFile']  , 'r' )    
    outputTEXT = open( parameters['OutputFile'] , 'w' )
    outputDiscardTEXT = open( Output_discard , 'w' )
    temp_recording_list = [  ]
    determination = None
    while True:
        try:
            line = next(inputTEXT).split('\n')[0].split('\r')[0]
            count_line = count_line + 1    
            ####
            if count_line % 2 == 1:
                temp_recording_list.append(    line    )
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
                if beginning == parameters['Marker']:
                    determination = 'yes'
                else:
                    determination = 'no'
                ##
                temp_recording_list.append(    line    )
                ## checking ##
                if len( temp_recording_list ) == 2:
                    pass
                else:
                    print    'Error ==> collecting: ' , temp_recording_list
                ######## Writing ########
                if determination == 'yes':
                    for written_line in temp_recording_list:
                        outputTEXT.write(  written_line + '\n'  )
                        count_available_line = count_available_line + 1    
                    count_available_reads = count_available_reads + 1    
                else:
                    for written_line in temp_recording_list:
                        outputDiscardTEXT.write(  written_line + '\n'  )
                        count_unavailable_line = count_unavailable_line + 1    
                    count_unavailable_reads = count_unavailable_reads + 1    
                ## deleting ##
                temp_recording_list = [  ]
                determination = None
        except StopIteration:
            break
    inputTEXT.close()    
    outputTEXT.close()
    outputDiscardTEXT.close()
    print    '\tChecking ==> Lines: ' , count_line , '    ==> ' , count_available_line + count_unavailable_line == count_line
    print    '\tChecking ==> Reads: ' , count_reads , '    ==> ' , count_reads , '* 2 =' , count_reads*2 , '    ==> ' , count_line == count_reads*2 , count_reads == sum( code_dict.values() ) , count_available_reads + count_unavailable_reads == count_reads
    print    '\tChecking ==> Available lines: ' , count_available_line
    print    '\tChecking ==> Available reads: ' , count_available_reads , '    ==> ' , count_available_reads , '* 2 =' , count_available_reads*2 , '    ==> ' , count_available_line == count_available_reads*2
    print    '\tChecking ==> Unavailable lines: ' , count_unavailable_line
    print    '\tChecking ==> Unavailable reads: ' , count_unavailable_reads , '    ==> ' , count_unavailable_reads , '* 2 =' , count_unavailable_reads*2 , '    ==> ' , count_unavailable_line == count_unavailable_reads*2
    ########
    Output_LOG = parameters['OutputFile'] + '.stats' + str(parameters['Number'])
    print    'Write down LOG: ' , Output_LOG
    outputTEXT = open( Output_LOG , 'w' )
    outputTEXT.write( 'Marker' +'\t'+ 'Count' +'\t'+ 'Percentage' +'\n' )
    for marker_count in sorted( code_dict.items() , key=lambda X:X[1] , reverse=True ):
        marker = marker_count[0]
        count  = str( marker_count[1] )
        percentage = str( float(count) / float(count_reads) )
        outputTEXT.write( marker +'\t'+ str(count) +'\t'+ str(percentage) +'\n' )
    outputTEXT.close()
    print    'Filtered ......'




if __name__ == '__main__':
    ####print    '\n'

    usage='usage: python R2_DNA_reads_filtering.py [options]'
    parser = OptionParser( usage = '%prog  -i /path/inputfile  -o /path/outputfile  -m nucleotide_marker  -n number_of_nucleotide' )
    parser.add_option('-i', '--input',  dest='input_file',                        type='string', help='"/path/filename"  Input file of fastq.')
    parser.add_option('-o', '--output', dest='output_file',                       type='string', help='"/path/filename"  Output file of filtered fastq.')
    parser.add_option('-m', '--marker', dest='nucleotide_marker_for_filter',      type='string', help='Nucleotides at the beginning of the read as the marker for filter.')
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
    if options.nucleotide_marker_for_filter == None:
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
    parameters[ 'Marker' ]     = options.nucleotide_marker_for_filter
    parameters[ 'Number' ]     = options.number_of_nucleotide_as_marker
    

    Fasta_filter( parameters )

    
    ####print    '\nThe end.\n'










