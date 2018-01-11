#June 2017 Stefano Bosisio
#Script to perform householder tridiagonalization on a matrix A

import os
import sys
from math import sqrt
import numpy as np



def norm(vect):
    val = 0.0
    for x in vect:
        val+=(x)**2

    return sqrt(val)

def sdiag(matrix):
    #extract sup and sub diagonal
    n = matrix.shape[0]
    supdiag = []  #super-diag
    subdiag = []  #sub-diag
    diag = [] # principal-diag
    for i in range(1,n):
        #matrix[row][col]
        row_sub = i
        col_sub = i-1
        row_sup = row_sub - 1
        col_sup = i
        sub = matrix[row_sub][col_sub]
        subdiag.append(sub)
        sup = matrix[row_sup][col_sup]
        supdiag.append(sup)
        princ = matrix[i][i]
        diag.append(princ)

    return subdiag,diag

def pythag(a,b):

    absa = abs(a)
    absb = abs(b)
    if (absa > absb):
        return absa*sqrt(1.0 + (absb/absa)*(absb/absa))
    else:
        return absb*sqrt(1.0 + (absa/absb)*(absa/absb))


def SIGN2(a,b):

    if b>=0.0:
        return -abs(a)
    else:
        return abs(a)
#take as an input a matix
#the matrix must be symmetric
#A =  np.array([[7,2,3,-1],[2,8,5,1],[3,5,12,9],[-1,1,9,7]],dtype=np.float64)
#A =  np.array([[4,1,-2,2],[1,2,0,1],[-2,0,3,-2],[2,1,-2,-1]],dtype=np.float64)
A =  np.array([[17,24,1,8,15],[23,5,7,14,16],[4,6,13,20,22],[10,12,19,21,3],[11,18,25,2,9]],dtype=np.float64)

AT = A.transpose()

#if np.allclose(A,AT):
#    print("The matrix is symmetric")
#else:
#    print("Error, the input matrix is not symmetric")
#    sys.exit(-1)

#matrix dimension
n = A.shape[0]

for k in range(n-2):
    #take the k-th column and delete the 0:k-th element
    u =(A[k][k+1:n]) #x vector, column and rows are the same since we have a symmetric matrix
#here we took the row [0] of matrix A (or row k) then A.A1 transform the row into an array
    #and we are taking the k+1 to n value of that array
    #create a vector u by computing the norm of x
    x = np.sign(u[0])*norm(u)
    u[0] = u[0]+x
    #compute the value H form the dot product u*u
    H = 0.5*np.dot(u,u)
    #REMEMBER with numpy array you have A[row:][:,column]
    #extract the submatrix taking out the k-th row and column
    subA = A[k+1:][:,k+1:n]
    #define a v vector as : subA*u / H
    v = subA.dot(u)/H
    #define a number g = dot(u,v)/(2.0*H)  = uT*v/2H
    g = np.dot(u.transpose(),v)/(2.0*H)
    #finally  v = v - g*u
    v = v - g*u
    #compute the transformation A'  = A' - vT*u - uT*v
    A[k+1:][:,k+1:n] = subA - np.outer(v,u) - np.outer(u,v)

    #eventually set A[k][k+1] = -u[0]
    A[k][k+1] = -x
    A[k+1][k] = -x
    #now store the diagonal, subdiagonal
print(A)
sub,diag= sdiag(A)

print("Diagonal elements after tridiagonalization\n")
print(diag)
print("Subdiag elements\n")
print(sub)
#TQLI


d = diag #main diagonalelements
e = sub # sub diagonal elements
#renumber the sub-diag elem e :
for i in range(1,len(e)):
    print("Substituting e[%d] = e[%d], %.2f and %.2f" %(i-1,i,e[i-1],e[i]))
    e[i-1] = e[i]
e.append(0.0)
print(e)
print(len(e))
print(len(d))
#now start with the iterative process
for l in range(0,len(d)):
    iteration = 0
    #now look for a single smal subdiagonal element to split the matrix
    for m in range(l,len(d)-1):
        dd = abs(d[m]) + abs(d[m+1])
        if float(abs(e[m])+dd) == dd:
            break #once it's found break this cycle

    if (m==l):
        break
    if (iteration ==30):
        print("ERROR too many iteration\n")
        sys.exit(-1)
    iteration +=1

    g = (d[l+1] - d[0])/(2.0*e[l])
    r = pythag(g,1.0)
    g = d[m] - d[l] + e[l] / (g+SIGN2(r,g))
    s = 1.0
    c = 1.0
    p = 0.0
    for i in range(m-1,l, -1):
        f = s*e[i]
        b = c*e[i]
        r = pythag(f,g)
        e[i+1] = r

        if (r==0.0):
            d[i+1] = d[i+1] - p
            e[m] = 0.0


        s = f/r
        c = g/r
        g = d[i+1] - p
        r = (d[i] - g)* s + 2.0*c*b
        p = s*r
        d[i+1] = g+p
        g = c*r - b
        #eigenvectors

    d[l] = d[l] - p
    e[l] = g
    e[m] = 0.0
    #
print("Eigenvalues")
print(d)
