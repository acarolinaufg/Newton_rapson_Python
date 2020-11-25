from Newton_rapson import newton

Exmplo = newton()

Exmplo.setBars(1,1,1.04,0.00,0+0*1j,0+0*1j)
Exmplo.setBars(2,2,1.00,0.00,1130e6+600e6*1j,0+0*1j)
Exmplo.setBars(4,2,1.00,0.00,650e6+100e6*1j,0+0*1j)
Exmplo.setBars(5,2,1.00,0.00,1000e6+300*1j,0+0*1j)
Exmplo.setBars(3,3,1.02,0.00,0+0*1j,1750e6+0*1j)

Exmplo.printBars()
Exmplo.setSesp()
Exmplo.ligaçoes(1,2,impedancia=0.0016+0.00672j)
Exmplo.ligaçoes(1,5,impedancia=0.00128+0.00504j)
Exmplo.ligaçoes(2,3,impedancia=0.00128+0.00504j)
Exmplo.ligaçoes(3,4,impedancia=0.00363+0.01344j)
Exmplo.ligaçoes(4,5,impedancia=0.00256+0.00608j)
Exmplo.ligaçoes(5,3,impedancia=0.00216+0.0044j)

#Exmplo.printLigaçoes()

#Exmplo.Ybus()

#Exmplo.Sinjetada()

#Exmplo.Jacobianas([2,4,5],[2,3,4,5])

Exmplo.solveCircuit(erro=0.0001,listTensao=[2,4,5],listAng=[2,3,4,5])
Exmplo.fluxoS(printTensao=True, printCorrentes=True)
Exmplo.Perdas()