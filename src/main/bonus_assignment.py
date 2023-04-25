import numpy as np
np.set_printoptions(precision=5, suppress=True, linewidth=100)
A = np.array([[3, 1, 1], [1, 4, 1], [2, 3, 7]])
b = np.array([1, 3, 0])
tol = 1e-6
iteration = 50
sum1=sum2=0
def norm(tol, x0, x):
    return(max(abs(x-x0)))/(max(abs(x0))+tol)
def seidel(A,b,iteration, tol):
    n= len(b)
    x= np.zeros(n)
    c=1
    while (c<= iteration):
        sum1=sum2=0
        x0 = x.copy()
        for i in range(n):
            for j in range(i):
                sum1 += A[i][j]*x[j]
            for j in range(i+1, n):
                sum2 += A[i][j]*x0[j]
            x[i]= (1/A[i][i])*(-sum1-sum2+b[i]) 
            if (norm(tol, x0,x)< tol):
                return c
        c += 1    
                
print(seidel(A,b,iteration, tol)) 
print()
def jacobi(A, b, iteration, tol):
    n = len(b)
    x = np.zeros(n)
    c = 1
    while c <= iteration:
        x0 = np.zeros(n)
        for i in range(n):
            total = 0
            for j in range(n):

                if (j != i):
                    total += A[i][j] * x[j]
            x0[i] = (b[i] - total) / A[i][i]
        if norm(tol,x0,x)<tol:
            return c
        x=x0
        c += 1
    return c
print(jacobi(A, b, iteration, tol))
print()
def f(x):
    return (x**3)- (x**2)+ 2

def prime(x):
    return 3* (x**2)- 2*x
x = 0.5
tol = 1e-6
iteration = 100
def newton(x,tol,iteration):
    k=0
    while abs(f(x))>= tol and k <iteration:
        quotient= f(x)/prime(x)
        x -= quotient
        k+=1
    print(k+1)  
    return k
newton(x,tol,iteration)
print()
def div_dif(matrix:np.array):
    size = len(matrix)
    for i in range(2,size):
        for j in range(2,i+2):
            if j>= len(matrix[i]) or matrix[i][j]!=0:
                continue
            
            left: float= matrix[i][j-1]
            diagonal: float= matrix[i-1][j-1]
            numerator: float= left-diagonal
            denominator = matrix[i][0]-matrix[i-j+1][0]
            operation= numerator/denominator
            matrix[i][j]=operation
    return matrix      
def hermite_interpolation():
    x_points = [0,1,2]
    y_points = [1,2,4]
    slopes = [1.06, 1.23, 1.55]
    num_of_points= len(x_points)
    matrix = np.zeros((2*num_of_points,2*num_of_points))
    i=0
    for x in range(0, num_of_points*2, 2):
        matrix[x][0]= x_points[i]
        matrix[x+1][0]=x_points[i]
        i+=1
    i=0    
    for x in range(0, num_of_points*2,2):
        matrix[x][1]=y_points[i]
        matrix[x+1][1]=y_points[i]
        i+=1
    i=0    
    for x in range(1, num_of_points*2,2):
        matrix[x][2]=slopes[i]
        i+=1
    filled_matrix= div_dif(matrix)
    print(filled_matrix)
hermite_interpolation() 
print()
def function(t,y):
    return y - t**3
def euler():
    a=0
    b=3
    n=100
    w=0.5
    step=(b-a)/n
    for i in range(n):
        w= w + step/2 * (function(a,w)+function(a+step,w+(step* function(a,w))))
        a+= step
    return round(w,5)
print(euler())  
print()