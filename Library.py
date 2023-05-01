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
        
def QR_EIG(A, A_prev):
    epsilon = 0.01
    Q,R = QR_Factor(A)
    sum = 0
    for i in range(len(A)):
        sum += A[i][i] - A_prev[i][i]
    if sum > epsilon:
        return QR_EIG(Mat_Mul(R,Q), A)
    else:
        return A
    
def SV_Mat(matrix):
    epsilon = 0.000001
    K = Mat_Mul(Transpose(matrix), matrix)
    A = QR_EIG(K,[[0 for _ in range(len(K[0]))] for _ in range(len(K))])
    SV = [[0 for _ in range(len(A[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        if A[i][i] > epsilon:
            SV[i][i] = A[i][i]
    return SV