import pandas as pd
import math


def Make_DF(filename):
    processed = []
    df = pd.read_csv(filename)
    for i, row in df.iterrows():
        processed.append(row["code"])
        processed.append(row["name"])
        processed.append(row["description"])
    return processed   


def Mat_Mul(A,B):
    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result    
              

def Transpose(A):
    result = [[0 for _ in range(len(A))] for _ in range(len(A[0]))]
    for i in range(0,len(A)):
        for j in range(0,len(A[0])):
            result[j][i]=A[i][j]
    return result

def QR_Factor(A):
    R = [[0 for _ in range(len(A[0]))] for _ in range(len(A))]
    Q = [[A[j][i] for i in range(len(A[0]))] for j in range(len(A))]
    for j in range(len(Q)):
        for i in range(len(Q)):
            R[j][j] += (Q[i][j])**2
        R[j][j] = math.sqrt(R[j][j])
        for i in range(len(Q)):
            Q[i][j]=Q[i][j]/R[j][j]
        for k in range(j+1, len(Q)):
            for i in range(len(Q)):
                R[j][k] += Q[i][j]*Q[i][k]
            for i in range(len(Q)):
                Q[i][k] = Q[i][k]-Q[i][j]*R[j][k]
    return(Q,R)

def QR_Eig(A, A_prev):
    epsilon = 0.0001
    Q,R = QR_Factor(A)
    sum = 0
    for i in range(len(A)):
        sum += abs(A[i][i]) - abs(A_prev[i][i])
    if sum > epsilon:
        return QR_Eig(Mat_Mul(R,Q), A)
    else:
        return A  
    
def Sing_Mat(matrix):
    epsilon = 0.000001
    K = Mat_Mul(Transpose(matrix), matrix)
    A = QR_Eig(K,[[0 for _ in range(len(K[0]))] for _ in range(len(K))])
    count = 0
    for i in range(len(A)):
        if A[i][i] > epsilon:
            count += 1
    SV = [[0 for _ in range(count)] for _ in range(count)]
    for i in range(count):
        if A[i][i] > epsilon:
            SV[i][i] = math.sqrt(A[i][i])
    return SV    

def Diag_Inv(A):
    inv = [[0 for _ in range(len(A[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        inv[i][i]=1/(A[i][i])
    return inv

#Returns the inverse of R, calling Inv_Col n times
def Up_Inv(R):
    n = len(R)
    R_inv = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        #basis vector
        e_i =[[0] for _ in range(n)]
        e_i[i][0] = 1
        col = Inv_Col(R, e_i)
        for j in range(n):
            #back substitution with ith basis vector
            R_inv[j][i] = col[j][0]

    return R_inv

#Returns column k of the inverse of R using back substitution with basis vectors
def Inv_Col(R, e_k):
    n = len(R)
    col_k = [[0] for _ in range(n)]

    for i in range(n-1, -1, -1):
        tmp = e_k[i][0]
        for j in range(n-1, i, -1):
            tmp -= col_k[j][0]*R[i][j]

        col_k[i][0] = tmp/R[i][i]
    return col_k


def Eig_Vect(A, e, uk):
    epsilon = 0.001
    K = Mat_Mul(Transpose(A),A)
    for i in range(len(K)):
        K[i][i] = K[i][i] - e
    Q,R = QR_Factor(K)
    Q_Inv = Transpose(Q)
    R_Inv = Up_Inv(R, )
    Coeff = Mat_Mul(R_Inv, Q_Inv)
    u = Mat_Mul(Coeff,uk)
    diff = 0
    for i in range(len(u)):
        diff += abs(u[i][0]) - abs(uk[i][0])
    if diff > epsilon:
        return Eig_Vect(A, e, u)
    else:
        return u

def SVD(A):
    S = Sing_Mat(A)
    V = [[0 for _ in range(len(S))] for _ in range(len(A[0]))]
    U = [[0 for _ in range(len(S))] for _ in range(len(A))]
    for i in range(len(V[0])):
        temp = Eig_Vect(A,S[i][i], [[1],[0],[0]])
        for j in range(len(V)):
            V[j][i] = temp[j][0]
    for i in range(len(U[0])):
        temp = Mat_Mul(A,[[V[j][i]] for j in range(len(V))])
        for j in range(len(U)):
            U[j][i] = temp[j][0]
    V_T = Transpose(V)
    return U,S,V_T

A = [[1,0,2],
     [3,1,0]]

U,S,V = SVD(A)

print(U)
print()
print(S)
print()
print(V)