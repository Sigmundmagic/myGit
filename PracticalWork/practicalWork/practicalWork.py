from math import sqrt

def splittingPDBString(string): # превращаем строку из PDB файла в массив строк. 'one two Z' -> ['one','two','Z']
    tmp = ''# временная строка
    arr = []# массив с результатом
    for i in range(len(string)):# решение для множественных пробелов
        if string[i] != ' ' :
            tmp += string[i]
        elif string[i] == ' ' and len(tmp) > 0 :
            arr.append(tmp)
            tmp = ''
    if tmp != '':
        arr.append(tmp)
    return arr

class managerForCalculatingDistancesBetweenAtoms:# объект, который сортирует аминокислоты в списке при их вставке
    def __init__(self):
        self.arrAtomAndAcidInf = [] # массив для объектов atomAndAcidInf
        self.informationAboutPositionOfOutermostAtom = atomAndAcidInf('','','',-1,0.0,0.0,0.0) # атом N-конца
    def addInfoFromLine(self,line): # приходит строка из файла
        res = splittingPDBString(line) 
        if len(res) != 12:# информация об атоме всегда содержит 12 элементов.Провека
            return 
        self.addInfo( atomAndAcidInf(res[2],res[3],res[4],int(res[5]),float(res[6]),float(res[7]),float(res[8])) )
    def addInfo(self,itemToAdd): # добавляет объекты с информацией об атомах в self.arrAtomAndAcidInf 
        if len(self.arrAtomAndAcidInf) == 0 and itemToAdd.dataType == 'N':
            self.informationAboutPositionOfOutermostAtom = itemToAdd
            return
        elif len(self.arrAtomAndAcidInf) > 0 and itemToAdd.dataType == 'N':
            return
        if self.informationAboutPositionOfOutermostAtom.dataType == '':
            return
        resultCalculateDistance = self.calculateDistance(itemToAdd)
        if resultCalculateDistance == -1:
            return
        itemToAdd.distanceToOutermostAtom = resultCalculateDistance 
        if len(self.arrAtomAndAcidInf) == 0   and itemToAdd.dataType == 'CA':
            self.arrAtomAndAcidInf.append(itemToAdd)
            return
        elif len(self.arrAtomAndAcidInf) != 0 and itemToAdd.dataType == 'CA':
            if self.arrAtomAndAcidInf[len(self.arrAtomAndAcidInf) - 1].distanceToOutermostAtom <= itemToAdd.distanceToOutermostAtom:
                self.arrAtomAndAcidInf.append(itemToAdd)#если расстояние до N-конца (последнего элемента в листе) меньше, либо равно расстоянию до N-конца в объекте, который поступил в качестве входного параметра, то добавляем этот объект в конец листа
            else:# если расстояние до N-конца в последнем элементе больше чем в поступающем, то мы удваиваем последний элемент и заменяем, а получившийся предпоследний заменяем на поступивший
                self.arrAtomAndAcidInf.append(self.arrAtomAndAcidInf[len(self.arrAtomAndAcidInf) - 1])
                self.arrAtomAndAcidInf[len(self.arrAtomAndAcidInf) - 2] = itemToAdd
    def calculateDistance(self,itemToAdd):
        if self.informationAboutPositionOfOutermostAtom.dataType != 'N': # если в классе нет информации о N-конце, то возращаем -1
            return -1
        return sqrt( pow(itemToAdd.X - self.informationAboutPositionOfOutermostAtom.X,2) + pow(itemToAdd.Y - self.informationAboutPositionOfOutermostAtom.Y,2) + pow(itemToAdd.Z - self.informationAboutPositionOfOutermostAtom.Z,2) )

class atomAndAcidInf:# информация о атоме из строки
    def __init__(self,dataType,aminoAcidName,sequenceName,aminoAcidIndexInSequence,X,Y,Z):
        self.dataType                 = dataType # N или CA
        self.aminoAcidName            = aminoAcidName # имя аминокислоты, которой принадлежит этот атом
        self.sequenceName             = sequenceName # имя последовательности, которой принадлежит этот атом из аминокислоты
        self.aminoAcidIndexInSequence = aminoAcidIndexInSequence
        self.X                        = X
        self.Y                        = Y
        self.Z                        = Z
        self.distanceToOutermostAtom  = 0.0 # расстояние до N-конца

class listForAminoAcidSequence:
    def __init__(self):
        self.arr = []# массив для хранения объектов aminoAcidSequence
    def addStr(self,seqObj): # наращивание одной из совпавших по имени последовательностей 
        for item in self.arr:
            if item.name == seqObj.name :
                item.sequence += seqObj.sequence
    def addItem(self,seqObj):# seqObj -> aminoAcidSequence 
        for item in self.arr:# перебор массива для хранения объектов aminoAcidSequence
            if item.name == seqObj.name:
                self.addStr(seqObj)
                return
        self.arr.append(seqObj)
    def addInfoOfAtom(self,line):# 
        arrFromLine = splittingPDBString(line)
        for item in self.arr:
            if item.name == arrFromLine[4]:
                if item.sequence.find(arrFromLine[3]) == -1:
                    return
                item.managerForCalculating.addInfoFromLine(line) #  наполняем информацией объект хранящий лист с atomAndAcidInf
                return
class aminoAcidSequence:# хранит последовательность, её имя, информацию об атомах 
    def __init__(self,string):
        self.name     = splittingPDBString(string)[2] # one two Z -> ['one','two','Z']
        self.sequence = ''
        self.managerForCalculating = managerForCalculatingDistancesBetweenAtoms()
        self.add(string)
    def add(self,string):# для наращивания последовательности
        if self.name == splittingPDBString(string)[2]:# фильтр для строк, которые не относятся к текущей последовательности
            for item in range(4,len(splittingPDBString(string))):# в 4 оказывается имя аминокислоты
                self.sequence += splittingPDBString(string)[item] + ' '
    def obtainSequenceOfAminoAcidsInOrderOfTheirSpatialDistanceFrom_N_terminalAminoAcid(self): #получить последовательность аминокислот в порядке их пространственного удаления от N-концевой аминокислоты
        res = ''
        for item in self.managerForCalculating.arrAtomAndAcidInf:# в листе из менеджера аминокислоты уже идут в порядке удаления 
            res += item.aminoAcidName + ' '
        return res



def whatIsThisLine(string): # определяем что это за строка 1 - SEQRES содержит последовательность; 2 - адрес крайнего атома ИЛИ адрес центрального атома
    if splittingPDBString(string)[0] == 'SEQRES' :
        return 1
    elif (splittingPDBString(string)[0] == 'ATOM' and splittingPDBString(string)[2] == 'N') or (splittingPDBString(string)[0] == 'ATOM' and splittingPDBString(string)[2] == 'CA'):
        return 2
    else:
        return -1

if __name__ == "__main__":
    listAminoAcid = listForAminoAcidSequence() # объект хранящий в себе последовательности и их имена. А также алгоритмы для наполнения информацией
    try:
        fo = open('result.txt','w')
    except Exception:
        print('err create result file!')
        sys.exit()
    try:
        f = open('1.pdb','r')
    except Exception:
        print('err open .pdb file!')
        sys.exit()
    try:
        line = f.readline().strip()
        if line == '':
            raise Exception('.pdb file is empty')
        if whatIsThisLine(line) == 1 :
            listAminoAcid.addItem( aminoAcidSequence(line) ) # создаёт или наращивает информацию о последовательности 
        while True:
            line = f.readline().strip()
            if len(line) == 0:
                break
            if whatIsThisLine(line) == 1 :
                listAminoAcid.addItem( aminoAcidSequence(line) ) # либо создаст новую последовательность, либо дополнит информацию
            elif whatIsThisLine(line) == 2:
                listAminoAcid.addInfoOfAtom(line)
        
        for item in listAminoAcid.arr:
            fo.write('>' + item.name + '\n');
            fo.write('' + item.sequence + '\n');
        fo.write('>Distance\n');
        for item in listAminoAcid.arr:
            fo.write('>' + item.name + '\n');
            fo.write('' + item.obtainSequenceOfAminoAcidsInOrderOfTheirSpatialDistanceFrom_N_terminalAminoAcid() + '\n');
    except Exception as e:
        print(e)
    f.close()
    fo.close()
    print('programm end')