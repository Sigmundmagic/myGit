import array
import re
import sys

class Pair:
    def __init__(self):
        self.one = ''
        self.two = ''

def complement(seq):
    basepairs = {"A": "T", "G": "C", "T": "A", "C": "G"}
    comp = ""
    for base in seq:
        comp += basepairs.get(base)
    return comp

def getDegeneratePrimer(filePathForDegeneratePrimer):
    degeneratePrimer = ''
    try:
        fileStreamForReadDegeneratePrimer = open(filePathForDegeneratePrimer,'r')
        degeneratePrimer = fileStreamForReadDegeneratePrimer.read()
        if degeneratePrimer == '':
            raise Exception('err read')
    except Exception:
        print('err read degeneratePrimer from file ' + filePathForDegeneratePrimer)
        sys.exit()
    return degeneratePrimer

# [ [A,B,C,D], [L,M,N] ]
def getPlotsCoincidingWithPrimer(arr,string):
    pattern = ''
    for i in arr:
        pattern += '['
        for item in i:
            pattern += item
        pattern += ']'
    pattern = pattern
    tmp = re.finditer(pattern,string)
    result = []
    for i in tmp:
        result.append([i.start(),i.end()])
    return result

def getArrayByPrimer(primer):
    result = []
    table = [
        ['A','A'],
        ['G','G'],
        ['T','T'],
        ['C','C'],
        ['R','A','G'],
        ['Y','C','T'],
        ['S','G','C'],
        ['W','A','T'],
        ['K','G','T'],
        ['M','A','C'],
        ['B','C','G','T'],
        ['D','A','G','T'],
        ['H','A','C','T'],
        ['V','A','C','G'],
        ['N','A','G','T','C']]
    for item in primer:
        flag = False
        for index in range(0,len(table)):
            if table[index][0] == item:
                flag = True
        if flag == False:
            raise Exception('primer validation error in func ' + getArrayByPrimer)
    for item in primer:
        for index in range(0,len(table)):
            if table[index][0] == item:
                tmpArr = []
                for i in range(1,len(table[index])):
                    tmpArr.append(table[index][i])
                result.append(tmpArr)
    return result

def resultStrings(firstPrimer,secondPrimer,sequenceName,sequence):
    result = []
    amplicons = []
    complementSequence = complement(sequence)
    occurrencesStringsInSequence = getPlotsCoincidingWithPrimer(getArrayByPrimer(firstPrimer),sequence)
    occurrencesStringsInComplementSequence = getPlotsCoincidingWithPrimer(getArrayByPrimer(secondPrimer),complementSequence)
    i = 0
    indexOccurrenceInSequence = 0
    while indexOccurrenceInSequence < len(occurrencesStringsInSequence):
        while i < len(occurrencesStringsInComplementSequence):
            if (occurrencesStringsInSequence[indexOccurrenceInSequence][0] - occurrencesStringsInComplementSequence[i][1] - 1) > 0:
                amplicons.append([occurrencesStringsInComplementSequence[i][1] - 1,sequence[occurrencesStringsInComplementSequence[i][1] - 1:occurrencesStringsInSequence[indexOccurrenceInSequence][0]]])
                i += 1
                break
            elif (occurrencesStringsInComplementSequence[i][0] - occurrencesStringsInSequence[indexOccurrenceInSequence][1] - 1) > 0:
                amplicons.append([occurrencesStringsInSequence[indexOccurrenceInSequence][1] - 1,sequence[occurrencesStringsInSequence[indexOccurrenceInSequence][1] - 1:occurrencesStringsInComplementSequence[i][0]]])
                i += 1
                break
            elif (occurrencesStringsInSequence[indexOccurrenceInSequence][0] - occurrencesStringsInComplementSequence[i][1] - 1) == 0 and (occurrencesStringsInComplementSequence[i][0] - occurrencesStringsInSequence[indexOccurrenceInSequence][1] - 1) == 0:
                i += 1
                break
            i += 1
        indexOccurrenceInSequence += 1
    result.append('>' + sequenceName)
    result.append(sequence)
    for item in amplicons:
        result.append(str(item[0]) + ' ' + item[1])
    return result
            

firstPrimer = getDegeneratePrimer(r'C:\Users\L\source\repos\Python\PythonApplicationTaskNumberThird\PythonApplicationTaskNumberThird\pr1.TXT')
secondPrimer = getDegeneratePrimer(r'C:\Users\L\source\repos\Python\PythonApplicationTaskNumberThird\PythonApplicationTaskNumberThird\pr2.TXT')
try:
    fo = open('result.txt','w')
except Exception:
    print('err create result file!')
    sys.exit()

try:
    f = open('1.fasta','r')
except Exception:
    print('err open fasta file!')
    sys.exit()
try:
    line = f.readline().strip()
    if line == '':
        raise Exception('fasta file is empty')
    pair = Pair()
    while True:
        if line.find('>') != -1:
            pair.one = line[line.find('>')+1:]
        else:
            pair.two = line
            tmpArray = resultStrings(firstPrimer,secondPrimer,pair.one,pair.two)
            for line_out in tmpArray:
                fo.write(line_out + '\n');
            pair = Pair()
        line = f.readline().strip()
        if len(line) == 0:
            break
    f.close()
    fo.close()
except Exception as e:
     print(e)