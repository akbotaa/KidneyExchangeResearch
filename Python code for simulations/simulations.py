import networkx as nx
import numpy as np
import random as rn
import matplotlib as mpl
import pylab
import csv

S=30 #for Simulations
s=0
T=4 #maximum number of periods

AvN=1000 #average number of pairs arriving in each period

a1=3
b1=2
p_O=0.374+0.066
p_A=0.357+0.063
p_B=0.085+0.015
p_AB=0.034+0.006

cost=0.2

NumbOfPairs=[0]*S

Matches1=[0]*S
Matches2=[0]*S
Matches3=[0]*S
SucMatches1=[0]*S
SucMatches2=[0]*S
SucMatches3=[0]*S



#randomly define blood type of the node based on BT dist'ns
def blood_type(pO, pA, pB, pAB):
    aux=rn.random()
    if aux<=pO:
        bt='O'
    elif aux<=pO+pA:
        bt='A'
    elif aux<=pO+pA+pB:
        bt='B'
    else:
        bt='AB'
    return bt 

while s<S:
    
    N=0
    t=0 #counter for periods
    
    G=nx.Graph()

    print "Generating data..."

    #iterate over time from 0 to T-1
    while t<T:

         numb=np.random.poisson(AvN)
     
         #iterate through pairs of period t
         for i in range(numb):

             dbt=blood_type(p_O, p_A, p_B, p_AB)
             pbt=blood_type(p_O, p_A, p_B, p_AB)
         
             #add p-d nodes and edge in-between
             G.add_node(N+2*i+1, type='d', btype=dbt, period=t)
             G.add_node(N+2*i+2, type='p', btype=pbt, period=t)
         
             #iterate through existing nodes and add edges where needed
             #randomly define weight, if edge is succesful during each period of its lifetime
             for n in G:
                 if G.node[n]['period']-t<5:
                     w=rn.betavariate(a1, b1)
             
                     if G.node[n]['type']=='p':
                 
                 
                         if dbt=='O':
                             e=(N+2*i+1,n)
                             G.add_edge(*e, wht=w, type=1, weight=1)
                             
                         elif dbt=='A':
                             if G.node[n]['btype'] in ['A', 'AB']:
                                 e=(N+2*i+1,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
     
                         elif dbt=='B':
                             if G.node[n]['btype'] in ['B', 'AB']:
                                 e=(N+2*i+1,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
    
                         else:
                             if G.node[n]['btype']=='AB':
                                 e=(N+2*i+1,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
    
                 
                     if G.node[n]['type']=='d':
                 
                         if pbt=='O':
                             if G.node[n]['btype']=='O':
                                 e=(N+2*i+2,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
        
                         elif dbt=='A':
                             if G.node[n]['btype'] in ['A', 'O']:
                                 e=(N+2*i+2,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
    
                         elif dbt=='B':
                             if G.node[n]['btype'] in ['B', 'O']:
                                 e=(N+2*i+2,n)
                                 G.add_edge(*e, wht=w, type=1, weight=1)
  
                         else:
                             e=(N+2*i+2,n)
                             G.add_edge(*e, wht=w, type=1, weight=1)
             
                     
                 
             ###redefine the edge between pair          
             if G.has_edge(N+2*i+1, N+2*i+2)==True:
                 G[N+2*i+1][N+2*i+2]['wht']=0
                 G[N+2*i+1][N+2*i+2]['weight']=0
                 G[N+2*i+1][N+2*i+2]['type']=0
             else:
                 G.add_edge(N+2*i+1, N+2*i+2, wht=0, weight=0, type=0)
         
         print(t)                
         t=t+1    
         
         N=N+2*numb    
  
    
    print "Running Algos..."

    #=================== Algorithms ===================
    #-- preparation finding true outcomes of matches --
    
    time=0

    while time<T:

         S_t=G.subgraph( [n for n,attrdict in G.node.items() if (0<=time-attrdict 
         ['period']<5)  ])
    
         #defining weights based on how old is the patient
         for n in S_t:
             if S_t.node[n]['type']=='p':
                 for m in S_t.neighbors(n):
                     S_t[n][m]['weight']=S_t[n][m]['wht']*(1-cost*(time-S_t.node[n]['period']))
                 
                     hah=np.random.binomial(1,S_t[n][m]['weight'])
                 
                     if time-S_t.node[n]['period']==0:
                         G[n][m]['s0']=hah
                     elif time-S_t.node[n]['period']==1:
                         G[n][m]['s1']=hah
                     elif time-S_t.node[n]['period']==2:
                         G[n][m]['s2']=hah
                     elif time-S_t.node[n]['period']==3:
                         G[n][m]['s3']=hah
                     else:
                         G[n][m]['s4']=hah
         print(time)             
         time=time+1


    #==================== Algorithm 2 =====================

    print "Running Algo 2"
    time=0
    G2=G.copy()
    counter2=0
    successful2=0
    while time < T:
     
         #algo 1 matching
         S_t=G2.subgraph( [n for n,attrdict in G.node.items() if (0<=time-attrdict 
         ['period']<5)  ])

         #defining weights based on how old is the patient
         for n in S_t:
             if S_t.node[n]['type']=='p':
                 for m in S_t.neighbors(n):
                     S_t[n][m]['weight']=S_t[n][m]['wht']*(1-cost*(time-S_t.node[n]['period']))
     
     
         M2=nx.max_weight_matching(S_t, maxcardinality=True)
     
     
         for n in S_t:
             if S_t.node[n]['type']=='p': #nechetniye
                 if n-M2[n]!=1: #if matched
                     G2.remove_node(n)
                     G2.remove_node(n-1)
   
                     counter2=counter2+1
                 
                     if time-S_t.node[n]['period']==0:
                         if G[n][M2[n]]['s0']==1:
                             successful2=successful2+1
                     elif time-S_t.node[n]['period']==1:
                         if G[n][M2[n]]['s1']==1:
                             successful2=successful2+1
                     elif time-S_t.node[n]['period']==2:
                         if G[n][M2[n]]['s2']==1:
                             successful2=successful2+1
                     elif time-S_t.node[n]['period']==3:
                         if G[n][M2[n]]['s3']==1:
                             successful2=successful2+1
                     else:
                         if G[n][M2[n]]['s4']==1:
                             successful2=successful2+1
                 
         print(time)
         time=time+1
     

    #==================== Algorithm 1 =====================

    print "Running Algo 1..."
    time=0
    G1=G.copy()
    counter1=0
    successful1=0

    while time < T:
     
         #algo 1 matching
         S_t1=G1.subgraph( [n for n,attrdict in G.node.items() if (0<=time-attrdict 
         ['period']<5)  ])
     
         M1=nx.max_weight_matching(S_t1, maxcardinality=True)
     
         for n in S_t1:
             if S_t1.node[n]['type']=='p': #nechetniye
                 if n-M1[n]!=1: #if matched
                     G1.remove_node(n)
                     G1.remove_node(n-1)
                 
                     counter1=counter1+1
                 
                     if time-S_t1.node[n]['period']==0:
                         if G[n][M1[n]]['s0']==1:
                             successful1=successful1+1
                     elif time-S_t1.node[n]['period']==1:
                         if G[n][M1[n]]['s1']==1:
                             successful1=successful1+1
                     elif time-S_t1.node[n]['period']==2:
                         if G[n][M1[n]]['s2']==1:
                             successful1=successful1+1
                     elif time-S_t1.node[n]['period']==3:
                         if G[n][M1[n]]['s3']==1:
                             successful1=successful1+1
                     else:
                         if G[n][M1[n]]['s4']==1:
                             successful1=successful1+1
     
         print(time)
         time=time+1


    #==================== Algorithm 3 =====================
    
    print "Running Algo 3..."
    time=0
    G3=G.copy()
    counter3=0
    successful3=0

    while time<T:
        Smean=G3.subgraph( [n for n,attrdict in G.node.items() if (0<=time-attrdict 
        ['period']<5)  ])
    
        #defining weights based on how old is the patient
        for n in Smean:
            if Smean.node[n]['type']=='p':
                for m in Smean.neighbors(n):
                    Smean[n][m]['weight']=Smean[n][m]['wht']*(1-cost*(time-Smean.node[n]['period']))
    
        md=a1-1/(a1+b1-2)
    
        for u,v in Smean.edges():
            if (Smean[u][v]['weight']<0.2 and Smean[u][v]['type']==1):
                Smean.remove_edge(u,v)
    
        M3=nx.max_weight_matching(Smean, maxcardinality=True)
    
        for n in Smean:
            if Smean.node[n]['type']=='p': #nechetniye
                if n-M3[n]!=1: #if matched
                    G3.remove_node(n)
                    G3.remove_node(n-1)
                
                    counter3=counter3+1
                
                    if time-Smean.node[n]['period']==0:
                        if G[n][M3[n]]['s0']==1:
                            successful3=successful3+1
                    elif time-Smean.node[n]['period']==1:
                        if G[n][M3[n]]['s1']==1:
                            successful3=successful3+1
                    elif time-Smean.node[n]['period']==2:
                        if G[n][M3[n]]['s2']==1:
                            successful3=successful3+1
                    elif time-Smean.node[n]['period']==3:
                        if G[n][M3[n]]['s3']==1:
                            successful3=successful3+1
                    else:
                        if G[n][M3[n]]['s4']==1:
                            successful3=successful3+1
        
        print (time)
    
        time=time+1
    
    NumbOfPairs[s]=len(G)/2
    
    Matches1[s]=counter1
    Matches2[s]=counter2
    Matches3[s]=counter3
    
    SucMatches1[s]=successful1
    SucMatches2[s]=successful2
    SucMatches3[s]=successful3
    
    print "S"
    print(s)
    s=s+1


fl = open('/Users/akbota/Documents/spring,2016/ECON605/Project! A+/Python code/results2.csv', 'w')

writer = csv.writer(fl)
writer.writerow(NumbOfPairs)
writer.writerow(Matches1)
writer.writerow(Matches2)
writer.writerow(Matches3)
writer.writerow(SucMatches1)
writer.writerow(SucMatches2)
writer.writerow(SucMatches3)

fl.close()