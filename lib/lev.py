#!/usr/bin/env python3
# Run mutect2 in tumor only mode
import sys
import os

def printMatrix(M):
        for row in M:
                print row
        print

def med(s, t):  
        k = len(s) + 1
        l = len(t) + 1

        M = [[0 for i in range(l)] for j in range(k)]
        MTrace = [["" for i in range(l)] for j in range(k)]

        M[0][0] = 0

#        print(k)
#        print(l)
#        print(M)
        for i in xrange(0, k):
                M[i][0] = i
                MTrace[i][0] = s[i-1]

        for j in xrange(0, l):
                M[0][j] = j
                MTrace[0][j] = t[j-1]

        MTrace[0][0] = "DONE"

        for i in xrange(1, k):
                for j in xrange(1, l):

                        sub = 1
                        sub_op = "sub"
                        if s[i-1] == t[j-1]:
                                # equality
                                sub = 0
                                sub_op = "eq"


                        # deletion
                        min_value = M[i-1][j] + 1
                        op = "del"
                        if min_value > M[i][j-1] + 1:
                                # insertion
                                min_value = M[i][j-1] + 1
                                op = "ins"
                        if min_value > M[i-1][j-1] + sub:
                                # substitution
                                min_value = M[i-1][j-1] + sub
                                op = sub_op


                        M[i][j] = min_value
                        MTrace[i][j] = op                        

#        print "final Matrix"
#        printMatrix(M)
#        printMatrix(MTrace)

############ MY PARTIAL SOLUTION


        def array_append(array,x,y):
            ops_string = MTrace[x][y]
            if ops_string == 'ins':
                print("->"+MTrace[0][y])
                array.append(("Insert",MTrace[0][y]))
            elif ops_string == 'sub':
                print(MTrace[x][0]+">"+MTrace[0][y])
                array.append(("Substitute",MTrace[x][0],MTrace[0][y]))
            elif ops_string == 'eq':
                array.append(("Equal",MTrace[x][0],MTrace[0][y]))
            elif ops_string == 'del':
                print(MTrace[x][0]+">-")
                array.append(("Delete",MTrace[x][0]))


        i = len(s)
        j = len(t)

        ops_array = []
        base = M[i][j]
        array_append(ops_array,i,j)


        while MTrace[i][j] != "DONE":
            base = M[i][j]
            local_min = min(M[i][j-1],M[i-1][j],M[i-1][j-1])
            if base == local_min:
                i = i - 1
                j = j - 1
                array_append(ops_array,i,j)
            elif M[i][j-1] < M[i-1][j]:
                j = j -1
                array_append(ops_array,i,j)
            elif M[i-1][j] < M[i][j-1]:
                i = i - 1
                array_append(ops_array,i,j)
            else:
                i = i - 1
                j = j - 1
                array_append(ops_array,i,j)

#        print ops_array
#########

        return M[k-1][l-1]      

med(sys.argv[1], sys.argv[2])
