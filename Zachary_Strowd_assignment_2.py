import argparse
#!/usr/bin/env python3
"""
=============================================================================
Title : Zachary_Strowd_assigment_2.py
Description : This script gets FNA data from a file and uses Needleman-Wunsch to find optimal alignment. Ouputs the score and aligned sequence into a file.
Author : Zachary Strowd
Date : 10/09/23
Version : 1.0
Usage : python3 Zachary_Strowd_assigment_2.py
Notes : Needs a input file to read sequences.
Python Version: 3.11.3
=============================================================================
"""

def getFNAdata(input_file):
    #Getting the data from the file
    dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            if '>' in line:
               l = line.strip()
               dict[l] = ''
               temp = l
            else:
                
                dict[temp] = dict[temp] + line.strip()
    return dict

def writeToFNA(dict, output_file, score):
    #Writing to the output file
    with open(output_file, 'w') as f:
        #parsing through the dictionary with the keys and values
        for key, value in dict.items():
            identLine = key + ";" +  " score=" + str(score) + "\n"
            f.write(identLine)
            sequenceLine = value
            #calling split80 to get the 80 character limit
            split80(sequenceLine, f)
    f.close()               

def split80(line, f):
    count = 0
    current = 0
    length = len(line)
    for char in line:
        current += 1
        count += 1
        if count == 80:
            f.write(char)
            f.write("\n")
            count = 0
        elif length == current:
            f.write(char)
            f.write("\n")
        else:
            f.write(char)

def readBlossumMatrix(file):
    blossumMatrix = {}

    with open(file) as f:
      matrix = f.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()

    for row in lines:
      scores = row.split()
      row_name = scores.pop(0)
      blossumMatrix[row_name] = {}

      for column_name in columns:
       blossumMatrix[row_name][column_name] = scores.pop(0)
   
    gap = int(blossumMatrix["-"]["A"])

    f.close()
    return blossumMatrix, gap

def needleman_wuncsh(blossumMatrix, sequences, gap):
    temp = []
    for i in sequences.keys():
        temp.append(i)
    sequence1 = sequences[temp[1]]
    sequence2 = sequences[temp[0]]
   
    gap_penalty = gap 
    gp = 0

    # Initialization of Score and Traceback matrix
    row = len(sequence1)
    col = len(sequence2)
    matrix = [[0 for i in range(col+2)] for j in range(row+2)]
    matrix[1][1] = 0

    traceback = [["-" for k in range(col+2)] for l in range(row+2)]
    traceback[1][1] = "C"

    # Adding necessary values to both matrix and traceback
    for i in range(2, row+2):
        gp += gap_penalty
        matrix[i][0] = sequence1[i-2]
        traceback[i][0] = sequence1[i-2]

        matrix[i][1] = gp
        traceback[i][1] = "U"

    gp = 0
    for j in range(2, col+2):
        gp += gap_penalty
        matrix[0][j] = sequence2[j-2]
        traceback[0][j] = sequence2[j-2]

        matrix[1][j] = gp
        traceback[1][j] ="L"
    
    #calculating the score matrix with Needleman-Wunsch
    for i in range(2, row+2):
        for j in range(2, col+2):
            top = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap
            diagonal = matrix[i-1][j-1] + int(blossumMatrix[matrix[i][0]][matrix[0][j]])
            matrix[i][j] = max(top, left, diagonal)
            #Filling in Traceback Matrix
            if matrix[i][j] == diagonal:
                traceback[i][j] = "D"
            elif matrix[i][j] == top:
                traceback[i][j] = "U"
            elif matrix[i][j] == left:
                traceback[i][j] = "L"

    finalScore = matrix[len(sequence1)+1][len(sequence2)+1]

    #getting strings
    newTraceBack = []
    for i in traceback[1:]:
        print(i)
        newTraceBack.append(i[1:])

    str1 = ""
    str2 =  ""
    i = len(sequence1)
    j = len(sequence2)
    traceback = newTraceBack
    #Calculating Output Alignment Sequences
    while i > 0 or j > 0:
        if traceback[i][j] == 'D':
            str1 += sequence1[i-1] 
            str2 += sequence2[j-1]
            i -= 1
            j -= 1
        elif traceback[i][j] == 'U':
            str1 += sequence1[i-1]
            str2 += '-'
            i -= 1
        else:
            str1 += '-' 
            str2 += sequence2[j-1]
            j -= 1
    # returning the reversed string
    return str1[::-1], str2[::-1], finalScore

def main():
    print("Assigment 2 :: R#11682551")

    #Setting up the Parsing for the input and output files
    parser = argparse.ArgumentParser(description='Parser')
    parser.add_argument('-i', type = str, help='Input file path')
    parser.add_argument('-o',type = str, help='Output file path')
    parser.add_argument('-s', type = str, help='Score matrix file path')
    args = parser.parse_args()
    
    input_file = args.i
    output_file = args.o
    blosum_matrix = args.s

    sequences = getFNAdata(input_file)
    blosumMat, gap = readBlossumMatrix(blosum_matrix)
    results = needleman_wuncsh(blosumMat, sequences, gap)

    # Getting keys and storing them in an array
    sequence1, sequence2, score = results
    temp = []
    for i in sequences:
        temp.append(i)
    
    sequences[temp[0]] = sequence2
    sequences[temp[1]] = sequence1
    print(sequences)
    
    #writing to file to complete
    writeToFNA(sequences, output_file, score)


if __name__ == "__main__":
    main()