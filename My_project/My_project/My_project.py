import random
import re


class Sequence:
    def __init__(self, intro_size=50, promoter="ACTG", terminator="ACTG"):
        self.promoter = promoter
        self.terminator = terminator
        self.first_chain = ""
        self.basepairs = {"A": "T", "G": "C", "T": "A", "C": "G"}
        self.second_chain = ""
        self.first_string(intro_size)
        self.complement()

    def Chargoff_rule(self, sequence):
        compl = ""
        for base in sequence:
            compl += self.basepairs.get(base)
        return compl


    def random_str(self, size):
        result = ""
        for i in range(size):
            choose = random.randint(1, 4)
            if choose == 1:
                result += "A"
            elif choose == 2:
                result += "T"
            elif choose == 3:
                result += "G"
            elif choose == 4:
                result += "C"
        while result.find("ACGT") != -1 or result.find("TGCA") != -1:
            result = result.replace("ACGT", "").replace("TGCA", "")
        return result

    def first_string(self, intro_size):
        string = "ACGT" + self.random_str(intro_size) + "ACGT" + self.random_str(intro_size) + "ACGT"
        self.first_chain = string

    def complement(self):
        comp = ""
        for base in self.first_chain:
            comp += self.basepairs.get(base)
        self.second_chain = comp


    def find_ind(self): # возвращает индексы всех ACGT
        f_ch = []
        first_chain = self.first_chain
        k = 0
        while k < len(first_chain):
            k = first_chain.find("ACGT", k)
            if k == -1:
                return f_ch
            else:
                f_ch.append(k)
                k += 1
        return f_ch


    def find_trp(self,my_arr,pattern_size = 4): # возвращает массив с индексами расположения транспозонов[[начало, конец],...]
        tmp_arr = [-1, -1]
        result = []
        i = 0
        while i < len(my_arr):
            if tmp_arr[0] == -1:
                tmp_arr[0] = my_arr[i]
            elif tmp_arr[1] == -1:
                tmp_arr[1] = my_arr[i] + pattern_size
                i -= 1
            elif tmp_arr[0] != -1 and tmp_arr[1] != -1:
                result.append(tmp_arr)
                tmp_arr = [-1, -1]
                i -= 1
            i += 1
        if tmp_arr[0] != -1 and tmp_arr[1] != -1:
            result.append(tmp_arr)
        return result


    def random_transposition(self, result_list ,left_to_right = False, exchange = False, both = False):
        trp_num = len(result_list) # количество транспозонов
        choose = random.randint(1,trp_num)
        random_trp = result_list[choose - 1] # массив из индексов начала и конца случайного транспозона
        breaks = ["ACGT----;TGCA----","ACGT----;----TGCA","ACG--T;T--GCA","AC-GT;TG-CA","A--CGT;TGC--A","----ACGT;TGCA----","----ACGT;----TGCA"]
        choose_1 = random.randint(1,7)
        random_break = breaks[choose_1-1] # строка со случайным разрывом
        choose_2 = random.randint(1,len(self.find_ind()))
        sites = self.find_ind() # все возможные сайты связывания
        site = sites[choose_2 - 1] # индекс одного случайного сайта связывания
        if (left_to_right == True and exchange == True) or (left_to_right == True and both == True) or (exchange == True and both == True):
            print("Выберите одну или ноль операций: left_to_right, exchange или both")
        elif left_to_right == False and exchange == False and both == False: # транспозон переместился без переворота или обмена
            if random_break == "ACGT----;TGCA----": # в данном случае транспозон переезжает без изменений
                self.first_chain = self.first_chain[:site + 2] + self.first_chain[random_trp[0]:random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACGT----;----TGCA": # достраивается
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + self.first_chain[random_trp[0]:random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACG--T;T--GCA": # достраивается
                self.first_chain = self.first_chain[:site + 2] + "ACGCGT" + self.first_chain[random_trp[0] + 4: random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "AC-GT;TG-CA": # перенос без изменений
                self.first_chain = self.first_chain[:site + 2] + self.first_chain[random_trp[0]:random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "A--CGT;TGC--A": # достраивается
                self.first_chain = self.first_chain[:site + 2] + "ACGCGT" + self.first_chain[random_trp[0] + 4: random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;TGCA----": # достраивается
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + self.first_chain[random_trp[0]:random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;----TGCA": # перенос без изменений
                self.first_chain = self.first_chain[:site + 2] + self.first_chain[random_trp[0]:random_trp[1]] + self.first_chain[site + 2:]
                self.complement()
            else:
                print("something is wrong with a break")
        elif left_to_right == True and exchange == False and both == False: # Перенос транспозона с переворотом
            t = self.first_chain[random_trp[0]:random_trp[1]]
            u = t[::-1] # транспозон в обратном порядке
            if random_break == "ACGT----;TGCA----":
                a = u[:len(u)-5] # изменившийся транспозон
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACGT----;----TGCA":
                a = u[:len(u) - 5]
                self.first_chain = self.first_chain[:site + 2] + "ACGTTGCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACG--T;T--GCA":
                a = u[:len(u) - 4]
                self.first_chain = self.first_chain[:site + 2] + "ACGAC" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "AC-GT;TG-CA":
                a = u[:len(u) - 3]
                self.first_chain = self.first_chain[:site + 2] + "AC" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "A--CGT;TGC--A":
                a = u[:len(u) - 2]
                self.first_chain = self.first_chain[:site + 2] + "ACG" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;TGCA----":
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + u + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;----TGCA":
                self.first_chain = self.first_chain[:site + 2] + u + self.first_chain[site + 2:]
                self.complement()
            else:
                print("something is wrong with a break")
        elif left_to_right == False and exchange == True and both == False: # вставка с обменом между цепями
            u = self.Chargoff_rule(self.first_chain[random_trp[0]:random_trp[1]])
            a = u[4:]
            if random_break == "ACGT----;TGCA----":
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACGT----;----TGCA":
                self.first_chain = self.first_chain[:site + 2] + "ACGTTGCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACG--T;T--GCA":
                self.first_chain = self.first_chain[:site + 2] + "ACGGCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "AC-GT;TG-CA":
                self.first_chain = self.first_chain[:site + 2] + "ACCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "A--CGT;TGC--A":
                self.first_chain = self.first_chain[:site + 2] + "ACGGCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;TGCA----":
                self.first_chain = self.first_chain[:site + 2] + "ACGTTGCA" + a + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;----TGCA":
                self.first_chain = self.first_chain[:site + 2] + u + self.first_chain[site + 2:]
                self.complement()
            else:
                print("something is wrong with a break")
        elif left_to_right == False and exchange == False and both == True: # Обмен с переворотом транспозона
            both_trp = self.Chargoff_rule(self.first_chain[random_trp[0]:random_trp[1]])[::-1]
            if random_break == "ACGT----;TGCA----":
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + both_trp[:len(both_trp)-5] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACGT----;----TGCA":
                self.first_chain = self.first_chain[:site + 2] + "ACGTACGT" + both_trp[:len(both_trp) - 5] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "ACG--T;T--GCA":
                self.first_chain = self.first_chain[:site + 2] + "ACGACG" + both_trp[:len(both_trp) - 4] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "AC-GT;TG-CA":
                self.first_chain = self.first_chain[:site + 2] + "AC" + both_trp[:len(both_trp) - 3] + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "A--CGT;TGC--A":
                self.first_chain = self.first_chain[:site + 2] + "ACG" + both_trp[:len(both_trp) - 2] + self.first_chain[site + 2:]
            elif random_break == "----ACGT;TGCA----":
                self.first_chain = self.first_chain[:site + 2] + "ACGT" + both_trp + self.first_chain[site + 2:]
                self.complement()
            elif random_break == "----ACGT;----TGCA":
                self.first_chain = self.first_chain[:site + 2] + both_trp + self.first_chain[site + 2:]
                self.complement()
            else:
                print("something is wrong with a break")





R1 = Sequence()
print(R1.__dict__)
R1.random_transposition(R1.find_trp(R1.find_ind()),)
print(R1.__dict__)

