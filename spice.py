import numpy as np
import os

##Component_Type | Node1 | Node2 | Value | Initial

### Function: Count number of nodes
t=0
result = dict()
def nodes_sources(lines):
    nodes_numbers = set()
    sources = 0
    for line in lines:
        x = line.split(" ")
        nodes_numbers.add(int(x[1][1]))
        nodes_numbers.add(int(x[2][1]))
        if 'Vsrc ' in line or 'I ' in line:
            sources+=1
    return max(nodes_numbers), sources

## File reading
fname = 'testcases/7.txt'
path, file = os.path.split(fname)
if not os.path.exists("output"):
    os.makedirs("output")
foutput = 'output/testcase_output_'+file
with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content] ##List of lines
step = float(content[0])
iterations = int(content[1])
content = content[2:]

##We gonna make the G, B, C, D Matrix
#Get number of nodes and Number of independent sources
np.set_printoptions(formatter={'float_kind':'{:f}'.format})
n, m = nodes_sources(content)
G = np.zeros((n,n))
B = np.zeros((n,m))
D = np.zeros((m,m))
Z = np.zeros((m+n, 1))
X = np.zeros((m+n, 1))
index_z = n
index_z_2 = 0
VCount = -1
caps = list()
inductors = list()
## Building G, B, D Matrices
for line in content:
    x = line.split(" ")
    node1 = int(x[1][1])
    node2 = int(x[2][1])
    value = float(x[3])
    init = int(x[4])

    ## Resistance case
    if x[0] =='R':
        if node2==0:
            G[node1-1, node1-1] += 1 / value
        else:
            G[node1-1, node2-1] += (-1/value)
            G[node2-1, node1-1] += (-1/value)
            G[node1-1,node1-1] += 1/value
            G[node2-1,node2-1] += 1/value

    ## Capacitance case
    if x[0] == 'C':
        if node2==0:
            caps.append((value, node1 - 1, 0))
            G[node1-1, node1-1] += value/step
            Z[node1-1, 0] += (value / step)*init
        else:
            caps.append((value, node1 - 1, node2 - 1))
            G[node1-1, node2-1] += (-1*(value/step))
            G[node2-1, node1-1] += (-1*(value/step))
            G[node1-1,node1-1] += value/step
            G[node2-1,node2-1] += value/step
            ##Effect on Z
            Z[node1-1,0] += (value/step)*init
            Z[node2-1,0] += -1*(value/step)*init

    ##Independent Voltage sources case
    if x[0]=="Vsrc":
        VCount+=1
        if node2==0:
            B[node1-1, VCount] += 1
        else:
            B[node1-1, VCount] += 1
            B[node2-1, VCount] += -1

        ##Number of rows in Z = m+n, so index 0 -> m+n-1
        Z[index_z,0] = value
        index_z+=1

    ##Current sources
    if x[0]=="Isrc":
        if node2==0:
            Z[node1-1, 0] += value
        else:
            Z[node1-1,0] += value
            Z[node2-1,0] += -1*(value)

    ##Inductors case
    if x[0]=='I':
        VCount += 1
        if node2==0:
            B[node1-1, VCount] += 1
        else:
            B[node1-1, VCount] += 1
            B[node2-1, VCount] += -1
        Z[index_z,0] += -1*(value/step)*init
        inductors.append((value,index_z))
        D[index_z-n, index_z-n] += -1*(value/step)
        index_z+=1

C = B.transpose()

## G, B, C, D Are done. let's build the A Matrix
A = np.zeros((m+n, m+n))

## Fill nxn by G
A[0:n, 0:n] = G
A[n:,0:n] = C
A[0:n, n:] = B
A[n:,n:] = D

################# All Matrices are done ############################
count_iterations = 0
print A
a_inverse = np.linalg.inv(A)
X = np.matmul(a_inverse,Z) ##Here we get initial solution of X
Z_init = Z.copy()
while count_iterations < iterations:
    print Z_init
    t+= step
    result[t]=X
    Z = Z_init.copy()
    #From X we want to pass capacitors and inductors to numerical solution function
    ##Capacitors
    cap = 0
    while cap < len(caps):
        C, n1, n2 = caps[cap]
        if n2==0:
            Z[n1,0] += (C/step)*X[n1]
        else:
            v = (C/step)*(X[n1]-X[n2])
            Z[n1,0]  += v
            Z[n2,0] -= v
        cap+=1

    ##Inductors
    ind = 0
    while ind < len(inductors):
        L, index = inductors[ind]
        Z[index,0] += X[index,0]*-1*(L/step)
        ind+=1

    X = np.matmul(a_inverse, Z)
    count_iterations+=1

##Write results to file in required format
f = open(foutput,'w')
### Print data of nodes
node = 1
key = step
while node <= n:
    key = step
    f.write('V' + str(node)+'\n')
    while key <= iterations*step:
        f.write(str(key) +"  "+ str(result[key][node-1][0])+'\n')
        key+=step
    f.write('\n')
    node+=1

##print currents
current = n

##print I_Vsrc s
while current < (m+n)-len(inductors):
    key = step
    f.write("I_Vsrc" + str(current-(n))+ '\n')
    while key <= iterations * step:
        f.write(str(key) + "  " + str(result[key][current][0])+'\n')
        key += step
    f.write('\n')
    current+=1

##print inductors
while current < X.shape[0]:
    key = step
    inductors_count = 0
    f.write("I_L" + str(inductors_count)+'\n')
    while key <= iterations * step:
        f.write(str(key) + "  " + str(result[key][current][0])+'\n')
        key += step
    f.write('\n')
    current+=1
    inductors_count+=1