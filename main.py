# Dame
import copy
import time


# Returneaza indicii posibilelor mutari
def indici_diagonale(poz, jucator):
    l = []
    # Daca e alb si nu e rege, merge doar in jos
    if jucator == 'a':
        if poz % 8 in range(1, 7):
            l.extend([poz + 7, poz + 9])
        if poz % 8 == 0:
            l.append(poz + 9)
        if poz % 8 == 7:
            l.append(poz + 7)
    if jucator == 'n':
        if poz % 8 in range(1, 7):
            l.extend([poz - 7, poz - 9])
        if poz % 8 == 0:
            l.append(poz - 7)
        if poz % 8 == 7:
            l.append(poz - 9)
    # Daca e rege, poate pe toate diagonalele
    if jucator == 'A' or jucator == 'N':
        if poz % 8 in range(1, 7):
            l.extend([poz - 7, poz - 9, poz + 7, poz + 9])
        if poz % 8 == 7:
            l.extend([poz - 9, poz + 7])
        if poz % 8 == 0:
            l.extend([poz - 7, poz + 9])
    # Iau doar pozitiile care se incadreaza in tabla
    l = [x for x in l if x >= 0 and x < 64]
    return l


def devine_rege(poz, jucator):
    if jucator == 'a' and poz in range(56, 64):
        return 'A'
    if jucator == 'n' and poz in range(0, 8):
        return 'N'
    return False


class Joc:
    NR_COLOANE = 8
    NR_LINII = 8
    SIMBOLURI_JUC = ['a', 'n']
    JMIN = ['a', 'A']
    JMAX = ['n', 'N']
    GOL = '.'

    def __init__(self, tabla=None):
        tabla_init = [
            '.', 'a', '.', 'a', '.', 'a', '.', 'a',
            'a', '.', 'a', '.', 'a', '.', 'a', '.',
            '.', 'a', '.', 'a', '.', 'a', '.', 'a',
            '.', '.', '.', '.', '.', '.', '.', '.',
            '.', '.', '.', '.', '.', '.', '.', '.',
            'n', '.', 'n', '.', 'n', '.', 'n', '.',
            '.', 'n', '.', 'n', '.', 'n', '.', 'n',
            'n', '.', 'n', '.', 'n', '.', 'n', '.',
        ]
        self.matr = tabla if tabla else tabla_init

    def final(self):
        nr_albe = 0
        nr_negre = 0
        for i in range(len(self.matr)):
            if self.matr[i] in Joc.JMIN:
                nr_albe += 1
            if self.matr[i] in Joc.JMAX:
                nr_negre += 1
        miscari_1, miscari_2 = self.mutari(Joc.JMIN)
        miscari_11, miscari_22 = self.mutari(Joc.JMAX)
        if not miscari_1 and not miscari_2 and not miscari_11 and not miscari_22:
            return False
        if nr_albe and not nr_negre:
            return Joc.JMIN
        if nr_negre and not nr_albe:
            return Joc.JMAX
        if nr_albe == 0 and nr_negre == 0:
            return 'remiza'
        return False

    # Returneaza cate piese sunt pe ultimul rand, pentru euristica
    def backrow(self, jucator):
        if jucator in ['n', 'N']:
            return [i for i in range(56, 64)]
        if jucator in ['a', 'N']:
            return [i for i in range(0, 8)]

    def mutari(self, jucator):
        l_mutari_1 = []
        l_mutari_2 = []
        for poz in range(len(self.matr)):
            if self.matr[poz].lower() in jucator:
                miscari_un_pas, miscari_doi_pasi = self.miscari_piesa_aleasa(self.matr[poz], self.matr, poz)
                i = 0
                # Cat timp pot sa fac miscari de doi pasi, gasesc unele noi si le sterg pe cele pe care le aveam(saritura multipla)
                while i < len(miscari_doi_pasi):
                    a, miscari_noi = self.miscari_piesa_aleasa(self.matr[poz], miscari_doi_pasi[i][0],
                                                               miscari_doi_pasi[i][1])
                    if (miscari_noi):
                        miscari_doi_pasi.pop(i)
                        i = i - 1
                        for (x, y) in miscari_noi:
                            miscari_doi_pasi.append((x, y))
                    i += 1

                if miscari_doi_pasi:
                    for (miscare, dif) in miscari_doi_pasi:
                        if (isinstance(miscare, list)):
                            l_mutari_2.append(Joc(miscare))
                else:
                    for miscare in miscari_un_pas:
                        l_mutari_1.append(Joc(miscare))
        return l_mutari_1, l_mutari_2

    def simbol_jucator(self, miscare, poz):
        return miscare[poz]

    # Returneaza miscarile posibile pentru o piesa aleasa(nu include hopuri multiple)
    def miscari_piesa_aleasa(self, jucator, tabla, poz):
        miscari_un_pas = []
        miscari_doi_pasi = []
        juc_opus = ['a', 'A'] if jucator.lower() == 'n' else ['n', 'N']
        posibile_diag = indici_diagonale(poz, jucator)
        for posibila in posibile_diag:
            if tabla[posibila] == Joc.GOL:
                miscare = copy.deepcopy(tabla)
                miscare[posibila] = devine_rege(posibila, jucator) if devine_rege(posibila, jucator) else jucator
                miscare[poz] = Joc.GOL
                miscari_un_pas.append(miscare)
            directie = posibila - poz
            if tabla[posibila] in juc_opus and (posibila + directie) in range(0, 64) and (posibila % 8 not in [0, 7]) \
                    and tabla[posibila + directie] == Joc.GOL:
                miscare = copy.deepcopy(tabla)
                miscare[posibila + directie] = devine_rege(posibila + directie, jucator) if devine_rege(
                    posibila + directie, jucator) \
                    else jucator
                miscare[posibila], miscare[poz] = Joc.GOL, Joc.GOL
                miscari_doi_pasi.append((miscare, posibila + directie))

        return miscari_un_pas, miscari_doi_pasi

    def miscare_valida(self, poz_init, poz_final):
        dif = poz_final - poz_init
        simbol = self.matr[poz_init]
        if simbol not in Joc.JMIN:
            return False
        indici_diag = indici_diagonale(poz_init, simbol)
        if dif in [-7, -9, 7, 9] and poz_final in indici_diag and self.matr[poz_final] == Joc.GOL:
            return True
        # Minim 7*2 maxim 9*2
        if dif in [-14, -18, 14, 18]:
            if self.matr[poz_init + dif // 2] in Joc.JMAX and self.matr[poz_final] == Joc.GOL:
                return True
        return False

    def fct_euristica(self):
        nr_jmin, nr_jmax = 0, 0
        backrow_jmin = self.backrow(Joc.JMIN)
        backrow_jmax = self.backrow(Joc.JMAX)
        for poz in range(len(self.matr)):

            # PRIMA EURISTICA
            #Creste scorul daca e rege, daca are piese pe backrow, daca e pe patratul din mijloc, sau daca e pe cele 2
            #Randuri din mijloc
            if self.matr[poz] in Joc.JMAX:
                if self.matr[poz].isupper():
                    nr_jmax += 7.75
                else:
                    nr_jmax += 5
                if poz in backrow_jmax:
                    nr_jmax += 4

                if poz in [27, 28, 35, 36]:
                    nr_jmax += 2

                if poz % 8 in [3, 4]:
                    nr_jmax += 0.5
            if self.matr[poz] in Joc.JMIN:
                if self.matr[poz].isupper():
                    nr_jmin += 7.75
                else:
                    nr_jmin += 5

                if poz in backrow_jmin:
                    nr_jmin += 4

                if poz in [27, 28, 35, 36]:
                    nr_jmin += 2

                if poz % 8 in [3, 4]:
                    nr_jmin += 0.5

            # A DOUA EURISTICA
            '''
            if self.matr[poz] in Joc.JMIN:
                nr_jmin+=1
            if self.matr[poz] in Joc.JMAX:
                nr_jmax+=1
            '''

        return nr_jmax - nr_jmin

    def estimeaza_scor(self, adancime):
        t_final = self.final()
        if t_final == Joc.JMAX:
            return (999 + adancime)
        elif t_final == Joc.JMIN:
            return (-999 - adancime)
        elif t_final == 'remiza':
            return 0
        else:
            return self.fct_euristica()

    def __str__(self):
        sir = ''
        for nr_col in range(self.NR_COLOANE):
            sir += str(nr_col) + ' '
        sir += '\n'

        for lin in range(self.NR_LINII):
            k = lin * self.NR_COLOANE
            sir += (" ".join([str(x) for x in self.matr[k: k + self.NR_COLOANE]]) + "\n")
        return sir


class Stare:
    ADANCIME_MAX = None

    def __init__(self, tabla_joc, j_curent, adancime, parinte=None, scor=None):
        self.tabla_joc = tabla_joc
        self.j_curent = j_curent

        # adancimea in arborele de stari
        self.adancime = adancime

        # scorul starii (daca e finala) sau al celei mai bune stari-fiice (pentru jucatorul curent)
        self.scor = scor

        # lista de mutari posibile din starea curenta
        self.mutari_posibile = []

        # cea mai buna mutare din lista de mutari posibile pentru jucatorul curent
        self.stare_aleasa = None

    def jucator_opus(self):
        if self.j_curent in Joc.JMIN:
            return Joc.JMAX
        else:
            return Joc.JMIN

    def mutari(self):
        l_stari_mutari1 = []
        l_stari_mutari2 = []
        l_mutari_1, l_mutari_2 = self.tabla_joc.mutari(self.j_curent)
        juc_opus = self.jucator_opus()
        if l_mutari_2:
            l_stari_mutari2 = [Stare(mutare, juc_opus, self.adancime - 1, parinte=self) for mutare in l_mutari_2]
        else:
            l_stari_mutari1 = [Stare(mutare, juc_opus, self.adancime - 1, parinte=self) for mutare in l_mutari_1]
        return l_stari_mutari1, l_stari_mutari2

    def __str__(self):
        sir = str(self.tabla_joc) + "(Juc curent: " + self.j_curent + ")\n"
        return sir


""" Algoritmul MinMax """


def min_max(stare):
    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
        return stare

    # calculez toate mutarile posibile din starea curenta
    l_stari_mutari1, l_stari_mutari2 = stare.mutari()
    if l_stari_mutari2:
        stare.mutari_posibile = l_stari_mutari2
    else:
        stare.mutari_posibile = l_stari_mutari1

    # aplic algoritmul minimax pe toate mutarile posibile (calculand astfel subarborii lor)
    mutari_scor = [min_max(mutare) for mutare in stare.mutari_posibile]

    if stare.j_curent == Joc.JMAX:
        # daca jucatorul e JMAX aleg starea-fiica cu scorul maxim
        stare.stare_aleasa = max(mutari_scor, key=lambda x: x.scor)
    else:
        # daca jucatorul e JMIN aleg starea-fiica cu scorul minim
        stare.stare_aleasa = min(mutari_scor, key=lambda x: x.scor)

    stare.scor = stare.stare_aleasa.scor
    return stare


def alpha_beta(alpha, beta, stare):
    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
        return stare

    if alpha >= beta:
        return stare  # este intr-un interval invalid deci nu o mai procesez

    l_stari_mutari1, l_stari_mutari2 = stare.mutari()
    if l_stari_mutari2:
        stare.mutari_posibile = l_stari_mutari2
    else:
        stare.mutari_posibile = l_stari_mutari1
    if stare.j_curent == Joc.JMAX:
        scor_curent = float('-inf')

        for mutare in stare.mutari_posibile:
            # calculeaza scorul
            stare_noua = alpha_beta(alpha, beta, mutare)

            if (scor_curent < stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor
            if (alpha < stare_noua.scor):
                alpha = stare_noua.scor
                if alpha >= beta:
                    break

    elif stare.j_curent == Joc.JMIN:
        scor_curent = float('inf')

        for mutare in stare.mutari_posibile:
            stare_noua = alpha_beta(alpha, beta, mutare)

            if (scor_curent > stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor

            if (beta > stare_noua.scor):
                beta = stare_noua.scor
                if alpha >= beta:
                    break

    stare.scor = stare.stare_aleasa.scor

    return stare


def citire(element):
    raspuns_valid = False
    while not raspuns_valid:
        if element == 'coloana':
            valoare = input("Coloana=")
            if valoare == "exit":
                return "exit"
            else:
                valoare = int(valoare)
                if valoare in range(0, 8):
                    raspuns_valid = True
                else:
                    print("Coloana nu este in intervalul corect")
        if element == 'linie':
            valoare = input("Linie=")
            if valoare == "exit":
                return "exit"
            else:
                valoare = int(valoare)
                if valoare in range(0, 8):
                    raspuns_valid = True
                else:
                    print("Linia nu este in intervalul corect")
    return valoare


def afis_daca_final(stare_curenta):
    if stare_curenta.tabla_joc.final():
        print("S-a terminat")
        return True
    else:
        return False


def main():
    # initializare algoritm
    raspuns_valid = False
    while not raspuns_valid:
        tip_algoritm = input("Algorimul folosit? (raspundeti cu 1 sau 2)\n 1.Minimax\n 2.Alpha-beta\n ")
        if tip_algoritm in ['1', '2']:
            raspuns_valid = True
        else:
            print("Nu ati ales o varianta corecta.")

    # initializare ADANCIME_MAX
    raspuns_valid = False
    while not raspuns_valid:
        n = input("Introdu numarul dificultatii:1.Usor 2.Mediu 3.Greu")
        if n.isdigit() and int(n) in range(1, 4):
            Stare.ADANCIME_MAX = int(n) * 2
            raspuns_valid = True
        else:
            print("Trebuie sa introduceti un numar intre 1,2 si 3!")

    # initializare jucatori
    [s1, s2] = Joc.SIMBOLURI_JUC.copy()  # lista de simboluri posibile
    raspuns_valid = False
    while not raspuns_valid:
        Joc.JMIN = str(input("Doriti sa jucati cu {} sau cu {}? ".format(s1, s2)))
        if (Joc.JMIN in Joc.SIMBOLURI_JUC):
            raspuns_valid = True
        else:
            print("Raspunsul trebuie sa fie {} sau {}.".format(s1, s2))
    Joc.JMAX = s1 if Joc.JMIN == s2 else s2

    # initializare tabla
    tabla_curenta = Joc()
    print("Tabla initiala")
    print(str(tabla_curenta))
    scor_calc=scor_jucator=0
    # creare stare initiala
    stare_curenta = Stare(tabla_curenta, Joc.SIMBOLURI_JUC[0], Stare.ADANCIME_MAX)
    mutari_Jucator = mutari_Calculator= 0
    timp_init = int(round(time.time() * 1000))
    while True:
        if (stare_curenta.j_curent == Joc.JMIN):
            # muta jucatorul
            raspuns_valid = False
            nr_pasi_efectuati = 0
            while not raspuns_valid:
                try:
                    miscare_1, miscare_2 = stare_curenta.tabla_joc.mutari(Joc.JMIN)
                    if miscare_2:
                        print("Obligat sa sari!")
                    if not miscare_2 and nr_pasi_efectuati:
                        raspuns_valid = True
                        break
                    coloana_init, linie_init = citire('coloana'), citire('linie')
                    if (coloana_init == "exit" and linie_init == "exit"):
                        print("Iesim din program.")
                        break
                    index_init = linie_init * Joc.NR_COLOANE + coloana_init
                    if stare_curenta.tabla_joc.matr[index_init] != stare_curenta.j_curent:
                        print("Nu este simbolul tau acolo!")
                    else:
                        print("Acum introdu unde ai vrea sa muti")
                        coloana_fin, linie_fin = citire('coloana'), citire('linie')
                        index_fin = linie_fin * Joc.NR_COLOANE + coloana_fin
                        if not stare_curenta.tabla_joc.miscare_valida(index_init, index_fin):
                            print("Miscare invalida! Reluam.")
                        else:
                            dif = index_fin - index_init
                            if dif not in [-14, -18, 14, 18] and miscare_2:
                                print("Trebuie sa sari!")
                                print(str(stare_curenta))
                            else:
                                stare_curenta.tabla_joc.matr[index_fin] = Joc.JMIN
                                stare_curenta.tabla_joc.matr[index_init] = Joc.GOL
                                if dif in [-14, -18, 14, 18]:
                                    stare_curenta.tabla_joc.matr[index_init + dif // 2] = Joc.GOL
                                    print(str(stare_curenta.tabla_joc))
                                    nr_pasi_efectuati += 1
                                else:
                                    raspuns_valid = True
                                    break

                except ValueError:
                    print("Coloana trebuie sa fie un numar intreg.")

            # dupa iesirea din while sigur am valida coloana
            # deci pot plasa simbolul pe "tabla de joc"

            # afisarea starii jocului in urma mutarii utilizatorului
            print("\nTabla dupa mutarea jucatorului")
            print(str(stare_curenta))
            scor_jucator = stare_curenta.scor
            # testez daca jocul a ajuns intr-o stare finala
            # si afisez un mesaj corespunzator in caz ca da
            if afis_daca_final(stare_curenta) or linie_init == "exit" and coloana_init == "exit":
                print("Final")
                break
            mutari_Jucator+= 1
            # S-a realizat o mutare. Schimb jucatorul cu cel opus

            stare_curenta.j_curent = stare_curenta.jucator_opus()


        # --------------------------------
        else:  # jucatorul e JMAX (calculatorul)
            # Mutare calculator

            # preiau timpul in milisecunde de dinainte de mutare
            t_inainte = int(round(time.time() * 1000))
            if tip_algoritm == '1':
                stare_actualizata = min_max(stare_curenta)


            else:  # tip_algoritm==2
                stare_actualizata = alpha_beta(-5000, 5000, stare_curenta)
            linie1 = -1
            coloana1 = -1

            stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc
            print("Tabla dupa mutarea calculatorului")
            print(str(stare_curenta))

            # preiau timpul in milisecunde de dupa mutare
            t_dupa = int(round(time.time() * 1000))
            print("Calculatorul a \"gandit\" timp de " + str(t_dupa - t_inainte) + " milisecunde.")
            scor_calc = stare_curenta.scor
            if afis_daca_final(stare_curenta):
                print("Final")
                break
            mutari_Calculator += 1
            # S-a realizat o mutare. Schimb jucatorul cu cel opus
            stare_curenta.j_curent = stare_curenta.jucator_opus()

    print(f"Calculatorul a facut {mutari_Calculator} miscari, Jucatorul a facut {mutari_Jucator} miscari")
    timp_final = int(round(time.time() * 1000))
    print(f"Timpul total de rulare este de:{timp_final - timp_init} ms")
    if scor_jucator and scor_calc:
        print(f"Scorul tau este:{scor_jucator}.\nScorul calculatorului este:{scor_calc}")
