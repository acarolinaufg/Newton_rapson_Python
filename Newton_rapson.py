##This code is the implementation for the newton-rapson problem

##Importing the modules
#cmth - Complex numbers, math - sin, cos..., numpy - matrix, matplotlib - visualization
import cmath as cmt 
import math as mt
import numpy as np
#import matplotlib.pyplot as plt

#Creating a class newton
class newton:
    def __init__(self):
        self.Sbase = 1000e6

        self.Sesp = dict()
        self.Sbarras = dict()
        self.Ligaçoes = dict()
        self.ybus = list()

        self.dados  = dict()
        self.tensaoPlot = dict()
        self.angPlot = dict()
        self.listAng = list()
        self.listTensao = list()


        self.nPV = int()
        self.nPQ = int()

        self.J1 = list()
        self.J4 = list()
        
        self.y = list()
        self.x = list()
        self.V = dict()
        self.I = dict()
        self.__fluxoS = dict()
        self.perdas = 0
        

    def setBars(self, barra,code, tensao, ang, carga, geraçao):
        """
        codes: 1- V e tetha; 2 - P e Q; 3 - P e V.

        ang - graus
        tensao - pu
        carga e geração - va

        To convert degree to rad
        To VA to pu

        """
        self.dados[barra] = {'code': code, 'tensao': tensao, 'ang':mt.radians(ang), 
                            'carga':(carga/self.Sbase), 'geraçao':(geraçao/self.Sbase)}
        
        self.tensaoPlot[barra] = [tensao]
        self.angPlot[barra] = [ang]
    
    def printBars(self):
        """
        Method to print the informations in the screen
        """
        print('\n\n=============Dados===========')
        print('Sbase = ', self.Sbase, ' VA')
        for i in self.dados:
            print(self.dados[i])
        print('=================================')


    def setSesp(self):
        """
        Method to calculate the specified power in each bar
        """

        for i in self.dados:
            if self.dados[i]['code']==2: ##PQ Bar
                self.Sesp[i] = {'Pesp': np.real(self.dados.get(i)['geraçao']-self.dados.get(i)['carga']),
                                'Qesp':float(
                                    np.imag(self.dados.get(i)['geraçao']-self.dados.get(i)['carga'])
                                    )}

            elif self.dados[i]['code']==3:##PV Bar
                self.Sesp[i] = {'Pesp': np.real(self.dados.get(i)['geraçao']-self.dados.get(i)['carga']),
                                'Qesp':float(
                                    np.imag(self.dados.get(i)['geraçao']-self.dados.get(i)['carga'])
                                    )}
        print('\n\n========================Sesp========================')
        print(self.Sesp, 'pu')
        print('====================================================')

    def ligaçoes(self, barra1,barra2, impedancia = None, admitancia = None):
        """
        Method to show the relation between the bars
        parameters in pu
        Here you can pass admtancia or impedancia
        """
        if impedancia is None:
            impedancia = 1/admitancia
        elif admitancia is None:
            admitancia = 1/impedancia
        else:
            return 'Error!!Please inform admitiancia or indutancia'
        self.Ligaçoes[(barra1,barra2)]= {'Impedancia': impedancia, 
                                        'Admitancia': admitancia}
    
    def printLigaçoes(self):
        print('\n\n=================Ligações=================')
        for i in self.Ligaçoes:
            print('ligação ',i,'\t',self.Ligaçoes[i])
        print('==============================================')

    def printYbus(self):
        print('\n\n===============Ybus=======================')
        for i in self.ybus: print(i)
        print('==============================================')

    def Ybus(self):
        #Here we are starting the matriz Ybus with zeros
        #The lenght of the matrix depends the quantity of bars
        self.ybus = np.ones((len(self.dados),len(self.dados)), dtype=complex)

        #this 'for' is for complete the ybus matrix
        for i in range(len(self.ybus)):
            lin = []
            for j in range(len(self.ybus)):
                if i==j:
                    lin.append(0)
                else:
                    if self.Ligaçoes.__contains__(tuple([i+1,j+1])):
                        lin.append(-self.Ligaçoes.get(tuple([i+1,j+1]))['Admitancia'])
                    elif self.Ligaçoes.__contains__(tuple([j+1,i+1])):
                        lin.append(-self.Ligaçoes.get(tuple([j+1,i+1]))['Admitancia'])
                    else:
                        lin.append(0)
            for j in range(len(self.ybus)):
                if i==j:
                    lin[j] = -1*sum(lin)

            self.ybus[i] = lin

        self.printYbus()

        ##This is a counter for we know how many PV and PQ bars there are in the problem
        for i in self.dados:
            if self.dados.get(i)['code']==2:
                self.nPQ +=1
            elif self.dados.get(i)['code']==3:
                self.nPV += 1
        
    def Sinjetada(self):
        #This function calculate the injected power in the system
        self.Sinj = dict()
        self.deltaPeQ = []
        self.ResiduoP = []
        self.ResiduoQ = []

        for i in self.dados:
            soma1 = []
            soma2 = []
            if self.dados[i]['code'] != 1:
                for j in self.dados:
                #In soma1 list we have the summation of active power
                    soma1.append(
                                abs(self.ybus[i-1][j-1])*
                                abs(self.dados.get(i)['tensao'])*
                                abs(self.dados.get(j)['tensao'])*  
                                mt.cos(np.angle(self.ybus[i-1][j-1]) 
                                - self.dados.get(i)['ang']
                                +self.dados.get(j)['ang'])
                    )

                #In soma20 list we have the summation of reactive power
                    soma2.append(
                                -abs(self.ybus[i-1][j-1])*
                                abs(self.dados.get(i)['tensao'])*
                                abs(self.dados.get(j)['tensao'])* 
                                mt.sin(np.angle(self.ybus[i-1][j-1]) 
                                - self.dados.get(i)['ang']
                                +self.dados.get(j)['ang'])* 1j)

                #Here we calculate the DeltaP and DeltaQ
                self.ResiduoP.append(np.real(
                    self.Sesp.get(i)['Pesp']-sum(soma1)))

                if self.dados[i]['code'] ==2:
                    self.ResiduoQ.append(np.imag(
                        self.Sesp.get(i)['Qesp']*1j - sum(soma2)))

        #Here we put the DeltaP and Delta Q in a list
        for i in range(len(self.ResiduoP)):
            self.deltaPeQ.append(self.ResiduoP[i])
        for i in range(len(self.ResiduoQ)):
            self.deltaPeQ.append(self.ResiduoQ[i])
        
        
        print('\n================DELTA P e DELTA Q================')
        for i in range(len(self.ResiduoP)):
            print('Delta P = ',self.ResiduoP[i])
        for i in range(len(self.ResiduoQ)):
            print('Delta Q - ',self.ResiduoQ[i])

    #To calculate the Jacobian sub-matrix J1
    def setJ1(self,listAng,nPQ,nPV):
        """
        listAng: list of angles to be calculate in the circuit (PV and PQ Bars)
        nPQ: number of PQ bars
        nPV: number of PV bars
        return: return the sub-matrix J1
        """
        self.J1 = np.ones((nPQ +nPV,nPQ+nPV))

        mainDiagonal = []
        outDiagonal = []

        for i in listAng:
            soma = []
            for j in range(1,len(self.dados)+ 1, 1):
                if i != j:
                    soma.append(
                        abs(self.ybus[i-1][j-1])*
                        abs(self.dados.get(i)['tensao'])*
                        abs(self.dados.get(j)['tensao'])*
                        cmt.sin(cmt.phase(self.ybus[i-1][j-1])-
                            self.dados.get(i)['ang']+
                            self.dados.get(j)['ang']
                    )
                    )
            mainDiagonal.append(sum(soma))

        for i in listAng:
            for j in listAng:
                if i != j:
                    outDiagonal.append(
                        -abs(self.ybus[i-1][j-1])*
                        abs(self.dados.get(i)['tensao'])*
                        abs(self.dados.get(j)['tensao'])*
                        cmt.sin(cmt.phase(self.ybus[i-1][j-1])-
                            self.dados.get(i)['ang']+
                            self.dados.get(j)['ang']
                    )
                )

        m = 0
        for i in range(len(listAng)):
            for j in range(len(listAng)):
                if i == j:
                   self.J1[i][j] = np.real(mainDiagonal[j])  
                else:
                    self.J1[i][j] = np.real(outDiagonal[m])
                    m+= 1
        #print('\nJ1 = \n',self.J1)

        return self.J1

    #def setJ2(self, listTensao,listAng, nPQ, nPV):
        
        #self.J2 = np.zeros((nPQ+nPV,nPQ))

        #mainDiagonal = []
        #outDiagonal = []

        #for i in listAng:
        #    soma = []
        #    a = 0
        #    for j in range(1,len(self.dados)+ 1, 1):
        #        if i != j:
        #            soma.append(
        #               abs(self.ybus[i-1][j-1])*
        #                abs(self.dados.get(j)['tensao'])*
        #                cmt.cos(cmt.phase(self.ybus[i-1][j-1])-
        #                    self.dados.get(i)['ang']+
        #                    self.dados.get(j)['ang']
        #            )
        #            )
        #    a = (2*abs(self.dados.get(i)['tensao'])*abs(self.ybus[i-1][i-1])*
        #            cmt.cos(cmt.phase(self.ybus[i-1][i-1])))

        #    mainDiagonal.append(a + sum(soma))

        #for i in listAng:
        #    for j in listTensao:
        #        if i != j:
        #            outDiagonal.append(
        #                abs(self.ybus[i-1][j-1])*
        #                abs(self.dados.get(i)['tensao'])*
        #                cmt.cos(cmt.phase(self.ybus[i-1][j-1])-
        #                    self.dados.get(i)['ang']+
        #                    self.dados.get(j)['ang']
        #            )
        #        )

        #m = 0
        #for i in range(nPQ+nPV):
        #    k = nPV
        #    for j in range(nPQ):
        #        if i < nPV:
        #           self.J2[i][j] = np.real(mainDiagonal[m])
        #           m+= 1  
        #        elif i >= nPV:
        #            if i - nPV == j:
        #                self.J2[i][j] = np.real(mainDiagonal[j+nPV])
        #                k +=1
        #            else:
        #                self.J2[i][j] = np.real(outDiagonal[m])
        #                m+= 1  

        #print('\nk = ', k, '\n')
        #print('\nJ2 = \n',self.J2)

        #return self.J2


    #def setJ3(self, listTensao,listAng, nPQ, nPV):
        
    #    self.J3 = np.zeros((nPQ,nPQ+nPV))

        #mainDiagonal = []
        #outDiagonal = []

        #for i in listAng:
            #soma = []
            #for j in range(1,len(self.dados)+ 1, 1):
                #if i != j:
                    #soma.append(
                        #abs(self.ybus[i-1][j-1])*
                        #abs(self.dados.get(i)['tensao'])*
                        #abs(self.dados.get(j)['tensao'])*
                        #cmt.cos(cmt.phase(self.ybus[i-1][j-1])-
                        #    self.dados.get(i)['ang']+
                        #    self.dados.get(j)['ang']
                    #)
                    #)

        #    mainDiagonal.append(sum(soma))

        #for i in listAng:
        #    for j in listTensao:
        #        if i != j:
        #            outDiagonal.append(
        #                -abs(self.ybus[i-1][j-1])*
        #                abs(self.dados.get(i)['tensao'])*
        #                abs(self.dados.get(j)['tensao'])*
        #                cmt.cos(cmt.phase(self.ybus[i-1][j-1])-
        #                    self.dados.get(i)['ang']+
        #                    self.dados.get(j)['ang']
        #            )
        #        )

        #m = 0
        #for i in range(nPQ):
            #for j in range(nPQ+nPV):
            #    if j < nPV:
            #       self.J3[i][j] = np.real(mainDiagonal[m])
            #       m+= 1  
            #    elif j >= nPV:
            #        if j - nPV == i:
            #            self.J3[i][j] = np.real(mainDiagonal[i+nPV])
            #        else:
            #            self.J3[i][j] = np.real(outDiagonal[m])
            #            m+= 1  

        #print('\nk = ', k, '\n')
    #    print('\nJ3 = \n',self.J3)

    #    return self.J3
         
    def setJ4(self, listTensao,listAng, nPQ, nPV):
        
        self.J4 = np.ones((nPQ,nPQ))

        mainDiagonal = []
        outDiagonal = []

        for i in listAng:
            soma = []
            a = 0
            for j in range(1,len(self.dados)+ 1, 1):
                if i != j:
                    soma.append(
                        abs(self.ybus[i-1][j-1])*
                        abs(self.dados.get(j)['tensao'])*
                        cmt.sin(cmt.phase(self.ybus[i-1][j-1])-
                            self.dados.get(i)['ang']+
                            self.dados.get(j)['ang']
                    )
                    )
            a = (2*abs(self.dados.get(i)['tensao'])*abs(self.ybus[i-1][i-1])*
                    cmt.sin(cmt.phase(self.ybus[i-1][i-1])))

            mainDiagonal.append(-a - sum(soma))

        for i in listAng:
            for j in listTensao:
                if i != j:
                    outDiagonal.append(
                        -abs(self.ybus[i-1][j-1])*
                        abs(self.dados.get(i)['tensao'])*
                        cmt.sin(cmt.phase(self.ybus[i-1][j-1])-
                            self.dados.get(i)['ang']+
                            self.dados.get(j)['ang']
                    )
                )

        m = 0
        for i in range(nPQ):
            for j in range(nPQ):
                if i == j:
                    self.J4[i][j] = np.real(mainDiagonal[j+nPV])
                else:
                    self.J4[i][j] = np.real(outDiagonal[m])
                    m+=1
                

        #print('\nk = ', k, '\n')
        #print('\nJ4 = \n',self.J4)

        return self.J4

    #In this method we dont need the jacobian matrix
    def Jacobianas(self, listTensao, listAng):
        """
        Método utilizado para calcular a matriz Jacobiana.
        :param listTensao: Lista de tensões a serem calculadas no circuito. (Barras PQ)
        :param listAng: Lista de ângulos a serem calculados no circuito. (Barras PQ e PV)
        Printa a matriz Jacobiana.
        """

        J1 = self.setJ1(listAng, self.nPQ, self.nPV)  # (nPQ  + nPV) X (nPQ + nPV)
        J4 = self.setJ4(listTensao, listAng, self.nPQ, self.nPV)  # (nPQ) X (nPQ)

        

        print('\n\n==================== SUBMATRIZES  JACOBIANAS: ===========================')
        print('\nJ1 = ')
        for i in J1: print(i)
        print('\nJ4 = ')
        for i in J4: print(i)
        print('========================================================================')

    def linerSystem(self):
        """
        This method is to solve the linear System
        LinearSystem will be used several times during the implementation so we need to
        restart de variables with zero to avoid problems.
        """
        self.y = []
        self.y = np.linalg.solve(self.J4,self.ResiduoQ)#V
        self.x = []
        self.x = np.linalg.solve(self.J1, self.ResiduoP)#theta

        itsOkay1 = np.allclose(np.dot(self.J4, self.y), self.ResiduoQ)
        print('\n\t DEU CERTO?', itsOkay1)
        itsOkay2 = np.allclose(np.dot(self.J1, self.x), self.ResiduoP)
        print('\n\t DEU CERTO?', itsOkay2)

        angulo = []
        tens = []
        for i in range(len(self.x)):
            angulo.append(self.x[i])
        for i in range(len(self.y)):
            tens.append(self.y[i])

        m=0
        for i in range (len(self.dados)):
            if self.dados.get(i+1)['code'] != 1:
                self.dados[i+1]['ang'] += float(np.real(angulo[m]))
                #self.angPlot[i+1].append(self.dados[i+1]['ang'])
                m+=1
        m=0
        for i in range (len(self.dados)):
            if self.dados.get(i+1)['code'] == 2: #PQ bar
                self.dados[i+1]['tensao'] += float(np.real(tens[m]))
                #self.tensaoPlot[i+1].append(self.dados[i+1]['tensao'])
                m+=1
    
    def NovaInjeçao(self):
        """
        Method to calculate the new Power Injection in the bars
        """
        self.Sbarras = dict()

        for i in self.dados:
            soma1 = []
            soma2 = []
            if self.dados[i]['code'] != 2:
                for j in self.dados:
                #In soma1 list we have the summation of active power
                    soma1.append(
                                abs(self.ybus[i-1][j-1])*
                                abs(self.dados.get(i)['tensao'])*
                                abs(self.dados.get(j)['tensao'])*  
                                mt.cos(np.angle(self.ybus[i-1][j-1]) 
                                - self.dados.get(i)['ang']
                                +self.dados.get(j)['ang'])
                    )

                #In soma20 list we have the summation of reactive power
                    soma2.append(
                                -abs(self.ybus[i-1][j-1])*
                                abs(self.dados.get(i)['tensao'])*
                                abs(self.dados.get(j)['tensao'])* 
                                mt.sin(np.angle(self.ybus[i-1][j-1]) 
                                - self.dados.get(i)['ang']
                                +self.dados.get(j)['ang'])* 1j)

                #Here we put the new powers in the dictionary
            if self.dados[i]['code'] == 1:
                self.Sbarras[i]={'P ':np.real(sum(soma1)), 'Q ': np.imag(sum(soma2))}
            elif self.dados[i]['code'] == 3:
                self.Sbarras[i] = {'P ': 0, 'Q ': np.imag(sum(soma2))}
            for i in self.Sbarras:
                print(self.Sbarras)

        for i in self.dados:
            if self.dados[i]['code'] == 1:
                self.dados[i]['geraçao'] = self.Sbarras.get(i)['P ']+self.Sbarras.get(i)['Q ']*1j
            elif self.dados[i]['code'] == 3:
                self.dados[i]['geraçao'] = np.real(self.Sbarras.get(i)['P '])+self.Sbarras.get(i)['Q ']*1j

    def solveCircuit(self, erro = None, iteraçoes = None, listTensao = None, listAng = None):
        """
        This method is to solve te circuit

        erro: parameter that will stop our method
        """
        self.listTensao = listTensao
        self.listAng = listAng
        self.count = 1

        self.Ybus()
        self.Sinjetada()
        self.Jacobianas(listTensao=self.listTensao,listAng=self.listAng)
        self.linerSystem()

        if iteraçoes is None and erro is not None:
            pEq = list(map(abs, self.deltaPeQ))
            #ver se o valor de P e Q são menores que o erro
            #quando não houver nenhum valor menor que o erro, as iterações param
            teste = list(map(lambda m: True if (m< erro) else False, pEq))
            stop = teste.count(False)
            while True:
                self.Sinjetada()
                self.Jacobianas(listTensao=self.listTensao,listAng=self.listAng)
                self.linerSystem()
                self.count +=1
                pEq = list(map(abs, self.deltaPeQ))
                teste = list(map(lambda m: True if (m< erro) else False, pEq))
                stop = teste.count(False)
                if stop ==0: #num de variaveis falsas
                    break
        elif iteraçoes is not None and erro is None:
             while self.count<iteraçoes:
                self.Sinjetada()
                self.Jacobianas(listTensao=self.listTensao,listAng=self.listAng)
                self.linerSystem()
                self.count +=1
                #pEq = list(map(abs, self.deltaPeQ))
                #teste = list(map(lambda m: True if (m< erro) else False, pEq))
                #stop = teste.count(False)
                
        self.NovaInjeçao()

        if iteraçoes is not None:
            print('\n=============Numero de iterações = ', self.count)
        elif erro is not None:
            print('\nConverviu para um erro de ',erro,' . ')
            print('Convergiu em ', self.count, 'iterações')

    def printTensao(self):
        print('\n==================Tensoes==========================')
        for i in self.V:
            print('Barra: \t', i, '\t Tensao: \t',self.V.get(i),'\t[pu]' )
    
    #here the Newton_rapson method was implemented, and we will calculate the Voltage
    #It can be showed in the screen or not, we just need pass the parameter print=True or False
    def Tensoes(self,print = None):

        self.V = dict()
        for i in self.dados:
            self.V[i] = cmt.rect( self.dados.get(i)['tensao'],
                                    self.dados.get(i)['ang'])
        if print:
            self.printTensao()

    def printCorrentes(self):
       
        print('============================ CORRENTES: =======================================')
        for i in self.I:
            print('Ligação: \t', i, '\tCorrente = \t', self.I.get(i), '\t[pu]')
        print('=============================================================================')
    
    def Correntes(self, print=None):
        """
        This method is to calculate the values of the currents in each line, 
        The calculus is made for all the bars, in the bars where there arent ligation 
        the result must be zero.

        the currents are calculated by the Voltage angles
        """
        self.I = dict()
        self.Tensoes(print=None)
        for i in self.dados:
            soma= []
            for j in self.dados:
                if i==j:
                    continue
                else:
                    self.I[(i,j)] = ((self.V.get(i)-self.V.get(j))*self.ybus[i-1][j-1])
                soma.append(((self.V.get(i)-self.V.get(j))*self.ybus[i-1][j-1]))
            self.I[(i,j)] = sum(soma)
        if print:
            self.printCorrentes()

    def fluxoS(self,printTensao=None, printCorrentes=None):
        """
        Calcular o fluxo de potencia em todas as ligações do sistema
        """
        self.__fluxoS = dict()
        self.Tensoes(print=printTensao)
        self.Correntes(print=printCorrentes)

        for i in self.I:
            a = i[0]
            self.__fluxoS[i] = -self.V.get(a)*np.conjugate(self.I.get(i))
        print('\n=============Fluxo de Potencia=====================')
        for i in self.__fluxoS:
            print('ligação:\t', i, '\tFluxo=', self.__fluxoS.get(i), '\t[pu]')
        print('==============================================================')

        for i in self.dados:
            if self.dados.get(i)['code'] != 2:
                self.dados[i]['geraçao'] = self.__fluxoS.get((i,i))
    
    def Perdas(self):
        """
        Calcula todas as perdas do circuito
        """
        self.perdas = 0
        perdas = []
        for i in self.__fluxoS:
            perdas.append(self.__fluxoS.get(i))
        self.perdas = sum(perdas)
        print('\n====================Perdas======================')
        print(self.perdas, '\t[pu]')




        

