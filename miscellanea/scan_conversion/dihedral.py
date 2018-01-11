import numpy
import math
import matplotlib.pyplot as plt


phi  = numpy.linspace(-3.14,3.14,100)
formula_uno = []
formula_due = []
formula_tre = []
formula_quattro= []

deconv=[]

for i in phi:
    calcolo_uno = 7.7843*(1+math.cos(i))
    calcolo_due = 6.3237*(1+math.cos(2*i))
    calcolo_tre = 11.7535*(1+math.cos(3*i+3.14)) 
    calcolo_quattro = 1.2379*(1+math.cos(4*i))
    sum  = calcolo_uno + calcolo_due + calcolo_tre + calcolo_quattro   

    formula_uno.append(calcolo_uno)
    formula_due.append(calcolo_due)
    formula_tre.append(calcolo_tre)
    formula_quattro.append(calcolo_quattro)
    deconv.append(sum)

line_uno=plt.plot(phi, formula_uno, 'r--',marker ='o')
line_due=plt.plot(phi, formula_due, 'b--',marker ='o')
line_tre=plt.plot(phi, formula_tre, 'g--', marker = 'o') 
line_quattro=plt.plot(phi, formula_quattro, 'y--', marker = 'o')
plt.ylabel("cosine ",fontsize = 18)   #maybe is better to give greater font size

plt.xlabel(r'$\phi$',fontsize= 18)

plt.savefig("cosine_modified.png")

plt.clf()


line_deco = plt.plot(phi, deconv, 'r--', marker = 'o')
plt.ylabel("Decon", fontsize = 18)
plt.xlabel(r'$\phi$', fontsize=18)
plt.savefig("deconv.png")
plt.clf()
