####### CVRP MIP using Gurobi ########

import sys
import numpy as np
from scipy.spatial.distance import cityblock
import datetime
import matplotlib.pyplot as plt
import random
import gurobipy as gp
import copy
import numpy.ma as ma
from itertools import permutations

def preprocessing(lines):
    courier=int(lines[0].rstrip('\n'))
    items=int(lines[1].rstrip('\n'))
    load=list(map(int, lines[2].rstrip('\n').split()))
    weight=list(map(int, lines[3].rstrip('\n').split()))
    xc=list(map(int, lines[4].rstrip('\n').split()))
    yc=list(map(int, lines[5].rstrip('\n').split()))
    xc = xc[-1:] + xc[:-1] 
    yc = yc[-1:] + yc[:-1] 

    load = sorted(load, reverse=True)
    for i in range(len(load)):
        if sum(load[0:i])>= sum(weight):
            load1 = load[0:i]
            load=load1
            courier = i
            break
    
    dist=np.zeros((len(xc),len(xc)))
    for i in range(len(xc)):
        for j in range(len(xc)):
            dist[i,j] = int(cityblock([xc[i],yc[i]],[xc[j],yc[j]]))
            dist=dist.astype(int)

    return courier,items,load,weight,dist

def load_assign(load, weight, dist, vehiclestart,items):
    #print("##### Dist assignment failed. Trying load assignment #####")
    j = 0
    while j < items:
        vehiclestart[np.argmax(weight)][0] = int(np.argmax(load) + 1)
        vehiclestart[np.argmax(weight)][1] = j  # affido pacco più grosso a load più grosso
        load[np.argmax(load)] -= max(weight)  #
        weight[np.argmax(weight)] = 0
        j += 1
    return vehiclestart

def dist_assign(l, w, d, last_item, items):
    i = 0
    load1 = copy.deepcopy(l)
    weight1 = copy.deepcopy(w)
    dist1 = copy.deepcopy(d)
    vehiclestart = np.zeros((items, 2))
    while i < items:
        m = np.argmax(dist1[items, :])
        if load1[np.argmax(load1)] == l[np.argmax(load1)]:
            if load1[np.argmax(load1)] >= weight1[m]:
                vehiclestart[m][0] = int(np.argmax(load1) + 1)
                vehiclestart[m][1] = i  # affido pacco più lontano a load più grosso
                last_item[np.argmax(load1)] = m
                load1[np.argmax(load1)] -= weight1[m]
                weight1[m] = 0
                dist1[:, m] = 0
        else:
            n = np.argmin(ma.masked_where(dist1[int(last_item[np.argmax(load1)]), 0:items] == 0, dist1[int(last_item[np.argmax(load1)]), 0:items]))
            print("Courier has already a package, nearest next package:", n)
            print("Weight of the nearest package:", weight1[n])
            if load1[np.argmax(load1)] >= weight1[n]:
                vehiclestart[n][0] = int(np.argmax(load1) + 1)
                vehiclestart[n][1] = i  # affido pacco più vicino al precedente
                last_item[np.argmax(load1)] = n
                load1[np.argmax(load1)] -= weight1[n]
                weight1[n] = 0
                dist1[:, n] = 0
            else:
                v = load_assign(l, w, d, np.zeros((items, 2)),items)
                return v
        i += 1
    return vehiclestart

def warm_reorder(K, vehiclestart, vehicle, items):
    for k in K:
        tmp = 0
        cnt = 0
        i = 0
        n = 0
        while i < vehiclestart.shape[0]:
            if vehiclestart[i, 0] == k:
                if tmp == 0 and vehiclestart[i, 1] == min(vehiclestart[vehiclestart[:, 0] == k, 1]) and vehiclestart[i, 1] != 1000:
                    vehicle[tmp, i + 1, k-1] = 1
                    tmp += 1
                    cnt = vehiclestart[i, 1]
                    vehiclestart[i, 1] = 1000
                    n=i
                for j in range(vehiclestart.shape[0]):
                    if vehiclestart[j, 0] == k and vehiclestart[j, 1] > cnt and vehiclestart[j, 1] != 1000 and vehiclestart[ j, 1] == min(vehiclestart[vehiclestart[:, 0] == k, 1]) and tmp >= 1 and i != j:
                        vehicle[n + 1, j + 1, k-1] = 1
                        cnt = vehiclestart[j, 1]
                        vehiclestart[j, 1] = 1000
                        n = j
                        i = -1
                        break
                    elif vehiclestart[j, 0] != k and j == vehiclestart.shape[0] - 1  and tmp != 0 and sum(vehicle[n + 1, :, k-1])==0: # and vehiclestart[j, 1] != 1000
                        vehicle[n + 1, 0, k-1] = 1
                        break
                    elif vehiclestart[i, 0] == k and j == vehiclestart.shape[0] - 1  and tmp != 0: # and vehiclestart[j, 1] != 1000
                        vehicle[n + 1, 0, k-1] = 1
                        break
            i += 1
        for s in range(1,items):
            if sum(vehicle[s,:,k-1])>=2:
                vehicle[s,0,k-1]=0
    return vehicle

def gurobi_model_2(N, V, K, A, q, C, dist, vehicle1):
    mdl = gp.Model("VRP")

    # Cost matrix
    c = { index: v for index, v in np.ndenumerate(dist) if index[0]!=index[1] }

    # Decision Variable
    x = mdl.addVars(A, vtype=gp.GRB.BINARY, name ="x")
    mdl.modelSense = gp.GRB.MINIMIZE

    # Constraints
    # 1) each node visited only once
    mdl.addConstrs(gp.quicksum(x[i,j,k] for k in K for i in V if i != j) == 1 for j in N)
    mdl.addConstrs(gp.quicksum(x[i,j,k] for k in K for j in V if j != i) == 1 for i in N)

    # 2) for each vehicle it can exit and enter the depot at max 1 time
    mdl.addConstrs(gp.quicksum(x[0,j,k] for j in N) == 1 for k in K)
    mdl.addConstrs(gp.quicksum(x[j,0,k] for j in N ) == 1 for k in K)

    # 3) num arcs in = num arcs out
    mdl.addConstrs((gp.quicksum(x[i,j,k] for i in V if i!=j) - gp.quicksum(x[j,i,k] for i in V if i!=j)) == 0 for j in V for k in K)

    # 4) load constraint
    mdl.addConstrs(gp.quicksum(x[i,j,k]*q[j] for i in V for j in N if j != i) <= C[k] for k in K)

    # Objective Function
    mdl.setObjective(gp.quicksum(x[i,j,k]*c[i,j] for i,j,k in A))

    # MIP start
    for i,j,k in A:
        x[i,j,k].Start = vehicle1[i,j,k]

    # Solve
    mdl.Params.TimeLimit = 60*5  #seconds
    mdl.Params.Heuristics = 0.1
    mdl.Params.MIPFocus = 2
    mdl.Params.LazyConstraints = 1

    return mdl, x

def gurobi_model_1(N, V, K, A, q, C, dist, vehicle1):
    mdl = gp.Model("VRP")

    # Cost matrix
    c = { index: v for index, v in np.ndenumerate(dist) if index[0]!=index[1] }

    # Decision Variables
    x = mdl.addVars(A, vtype=gp.GRB.BINARY, name ="x")
    u = mdl.addVars(N, vtype=gp.GRB.CONTINUOUS, name ="u")
    mdl.modelSense = gp.GRB.MINIMIZE

    # Constraints
    # 1) each node visited only once
    mdl.addConstrs(gp.quicksum(x[i,j,k] for k in K for i in V if i != j) == 1 for j in N)
    mdl.addConstrs(gp.quicksum(x[i,j,k] for k in K for j in V if j != i) == 1 for i in N)

    # 2) for each vehicle it can exit and enter the depot at max 1 time
    mdl.addConstrs(gp.quicksum(x[0,j,k] for j in N) == 1 for k in K)
    mdl.addConstrs(gp.quicksum(x[j,0,k] for j in N ) == 1 for k in K)

    # 3) num arcs in = num arcs out
    mdl.addConstrs((gp.quicksum(x[i,j,k] for i in V if i!=j) - gp.quicksum(x[j,i,k] for i in V if i!=j)) == 0 for j in V for k in K)

    # 4) load constraint
    mdl.addConstrs(gp.quicksum(x[i,j,k]*q[j] for i in V for j in N if j != i) <= C[k] for k in K)

    # 5) Miller-Tucker-Zemlin formulation
    mdl.addConstrs((x[i,j,k]==1)>>(u[i]+q[j]==u[j]) for i in N for j in N for k in K if i!=j)
    mdl.addConstrs(q[i] <= u[i] for i in N)
    mdl.addConstrs((x[i,j,k]==1)>>(u[i] <= C[k]) for k in K for i in N for j in V if i!=j )

    # Objective Function
    mdl.setObjective(gp.quicksum(x[i,j,k]*c[i,j] for i,j,k in A))

    # MIP start
    for i,j,k in A:
        x[i,j,k].Start = vehicle1[i,j,k]

    # Solve
    mdl.Params.TimeLimit = 60*5  #seconds
    mdl.Params.Heuristics = 0.1
    mdl.Params.MIPFocus = 2

    return mdl, x


def plot_sol(A, x, courier, xc, yc, choice, argv):
    active_arcs = [a for a in A if x[a].x > 0.99]
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                for i in range(courier)]
    for i, j, k in active_arcs:
        plt.plot([xc[i], xc[j]], [yc[i], yc[j]], c=color[k-1], zorder=0)
    plt.plot(xc[0], yc[0], c='r', marker='s')
    plt.scatter(xc[1:], yc[1:], c='b', s=2.5)
    if choice=="1":
        plt.savefig('../plots/MTZ/route_'+argv+'.png')
    if choice=="2":
        plt.savefig('../plots/GSEC/route_'+argv+'.png')
    return

def tours_writer(out, K, active_arcs):
    for k in K:
        l=[0]
        i=0
        for tuple1 in active_arcs:
            for tuple in active_arcs:
                if tuple[0]==i and tuple[2]==k and tuple[1] not in l:
                    i=tuple[1]
                    l.append(tuple[1])
        l.append(0)
        print("Courier "+str(k)+": ")
        out.write("Courier "+str(k)+": ")
        for i in range(len(l)-1):
            print(str(l[i])+"->", end="")
            out.write(str(l[i])+"->")
        print(0, end="")
        out.write(str(0))
        print("\n")
        out.write("\n")
    out.close()
    return


def main(argv):

    # Callback - use lazy constraints to eliminate sub-tours
    def subtourelim(model, where):
        if where == gp.GRB.Callback.MIPSOL:
            # make a list of edges selected in the solution
            vals = model.cbGetSolution(model._x)
            selected = gp.tuplelist((i,j) for i, j, k in model._x.keys()
                                    if vals[i, j, k] > 0.5)
            # find the shortest cycle in the selected edge list
            tour = subtour(selected)
            if len(tour) < len(V): 
                for k in K:
                    model.cbLazy(gp.quicksum(model._x[i, j, k]
                                            for i, j in permutations(tour, 2))
                                <= len(tour)-1)


    # Given a tuplelist of edges, find the shortest subtour not containing depot (0)
    def subtour(edges):
        unvisited = list(range(1, len(V)))
        cycle = range(len(V)+1)  # initial length has 1 more city
        while unvisited:
            thiscycle = []
            neighbors = unvisited
            while neighbors:
                current = neighbors[0]
                thiscycle.append(current)
                if current != 0:
                    unvisited.remove(current)
                neighbors = [j for i, j in edges.select(current, '*')
                            if j == 0 or j in unvisited]
            if 0 not in thiscycle and len(cycle) > len(thiscycle):
                cycle = thiscycle
        return cycle

    data_file= open('./MCP_Instances/'+argv)
    lines=[]
    for line in data_file:
      lines.append(line)
    data_file.close()

    courier,items,load,weight,dist = preprocessing(lines)

    # Sets and parameters

    N = [i+1 for i in range(items)]
    V = [0]+N
    K = [i+1 for i in range(courier)]
    A = [(i,j,k) for i in V for j in V for k in K if i!=j]
    q = {key:weight[key-1] for key in N}
    C = {key:load[key-1] for key in K}

    last_item = np.zeros(courier)
    vehiclestart = dist_assign(load, weight, dist, last_item,items)
    for i in range(vehiclestart.shape[0]):
        vehiclestart[i, 0] = int(vehiclestart[i, 0])
    vehicle = np.zeros((items+1, items+1, courier))
    vehicle = warm_reorder(K,vehiclestart,vehicle,items)
    vehicle_final = {(i, j, k): int(vehicle[i,j,k-1]) for i, j, k in A}

    choice = input('Select model (1=MTZ, 2=GSEC):')

    if choice=="1":
        m,x  = gurobi_model_1(N,V,K,A,q,C,dist,vehicle_final)
        m.Params.Threads = 2
        m.optimize()
        out = open('../out/MTZ/out'+argv[1:][3:]+'.txt', 'w')
    elif choice=="2":
        m,x  = gurobi_model_2(N,V,K,A,q,C,dist,vehicle_final)
        m._x = x
        m.Params.Threads = 2
        m.optimize(subtourelim)
        out = open('../out/GSEC/out'+argv[1:][3:]+'.txt', 'w')

    try:
        pass
    except AttributeError as error:
        print('Solution not found')

    print("Solution found")
    out.write("Solution found")
    out.write("\n")
    print('Total distance: %g' % m.ObjVal)
    out.write('Total distance: %g' % m.ObjVal)
    out.write("\n")

    active_arcs = [a for a in A if x[a].x > 0.99] 
    tours_writer(out, K, active_arcs)

    xc=list(map(int, lines[4].rstrip('\n').split()))
    yc=list(map(int, lines[5].rstrip('\n').split()))
    xc = xc[-1:] + xc[:-1] 
    yc = yc[-1:] + yc[:-1]

    try:
        plot_sol(A, x, courier, xc, yc, choice, argv)
    except:
        print("No solution, plotting MIP start...")
        plot_sol(A, vehicle_final, courier, xc, yc, choice, argv)

if __name__ == "__main__":
    main(sys.argv[1])