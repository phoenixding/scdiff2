#!/usr/bin/env python
# coding: utf-8

# In[1]: Author statement


#!/usr/bin/env python
# Author: Jun Ding
# Email: junding (at) cs (dot) cmu (dot) edu
# Date: June. 29th, 2020
# 
# This scdiff software suite is desinged to infer the clusters, trajectories, and regulatory
# networks underlying dynamic biological process (e.g., cell differntiation, disease progression)
# based on given time-series single-cell expression input data.  Please use "scdiff -h" for the detailed usage.
# 
# This software is freely avaible for academic uses. 
# For any commerical usage, please contact me at the email address above.
# All rights reserved.
# Please don NOT modify the above statement.

# In[2]: Import modules
## Log in needed modules 
    
import pdb,sys,os,random
import numpy as np
import scanpy as sc
import anndata
import argparse
import datetime
import functools
import math
from multiprocessing import Pool
import warnings
warnings.simplefilter("ignore")

from File import * 
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from scipy.stats import zscore
from scipy.stats import norm
from scipy.stats import binom
from scipy.stats import ranksums
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
from scipy.sparse import csr_matrix,csc_matrix,lil_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from imblearn.over_sampling import SMOTE

from functools import partial

import pkg_resources
from viz2 import *

# In[3]:Class cell

## class cell
class Cell:
    def __init__(self, Cell_ID, TimePoint, Expression,typeLabel):
        self.ID=Cell_ID
        self.T=TimePoint
        self.E=Expression
        self.Label=None # Label for clustering purpose
        self.typeLabel=typeLabel
        



# In[4]:Class cluster(Node)
## class Cluster (Node in the Tree Graph-Trajectory)
class Cluster:
    def __init__(self,cells,ID):
        self.cells=cells                                 # cells (data poitns) for the cluster
        self.ID=int(ID)                                  # ID (must be integer=> matches the cell labels from the pre-run results)
        self.P=None                                      # parent cluster
        self.C=[]                                        # child clusters
        self.EMatrix=self.__getEMatrix()               # get sparse (csc sparse) expression matrix for all the cells in the Cluster 
        [self.E,self.R]=self.__getAvgVarE()              # initial mean,variance expression
        [self.mT,self.rT]=self.__getAvgVarT()            # initial mean,variance Time
        self.T=self.mT                                   # Time 
        self.PR=0                                        # Prior probability of cluster
            
    
    #--------------------------------------------------------------
    ## public functions:
    
    # calculate the total probability of a cell belonging to given cluster (based on expression, and time)
    def getAssignProbability(self,cell,W,TW=0.5):
        # public interface to calculate probability of cell assignment 
        # W: Gene weight
        # TW: time weight
        #pdb.set_trace()
        print("assign cell %s"%(cell.ID))
        PG = self.__getCellProbability(cell,W) 
        PT=TW*norm.logpdf(cell.T,self.mT,self.rT)
        PR = math.log(self.PR)
        P=PG+PT+PR
        return P
    
    # get expressed TFs of the cluster
    def getExpressedTF(self,TFList,GL):
        print("get expressed TFs in node %s ..."%(self.ID))
        PTF=[]
        EXCUT=0.15   # at least 1-EXCUT non-zero => expressed
        indptr=self.EMatrix.indptr
        M=len(self.cells)
        
        for i in TFList:
            if i in GL:
                ix=GL.index(i)
                nonzero=indptr[ix+1]-indptr[ix]
                if nonzero*1.0/M>EXCUT:
                    PTF.append(i)
        return PTF
    #----------------------------------------------------------------
    ## private functions
    # calculate the probability of a cell belonging to given cluster  based on the guassian distribution (expression)
    def __getCellProbability(self,cell,W):
        # priviate function, can't be accessed outside, can only be called by other functions of the class
        # cell: a cell
        # mu,sm
        # return mean logpdf (log of the probability density function) => log joint density 
        mu=self.E
        sm=self.R
        N=len(mu)
        cx=cell.E.toarray()[0]
        #pdb.set_trace()
        P=[]
        for i in range(N):
            p=math.log(W[i])+norm.logpdf(cx[i],mu[i],sm[i])
            P.append(p)
        return np.mean(P)
    
    # get the mean,sigma for the time of cells in the cluster
    def __getAvgVarT(self):
        iT=[item.T for item in self.cells]
        pCount=0.01
        imu=round(np.mean(iT),2)
        ivar=round(np.var(iT)+pCount,2)
        return [imu,ivar]
        
    # get the sparse cell expression matrix
    # return a private expression matrix (sparse)
    def __getEMatrix(self):
        print("get sparse expression matrix for %s ..."%(self.ID))
        M=len(self.cells)
        N=self.cells[0].E.shape[1]
        mtx=lil_matrix((M,N)) # lil_matrix (sparse) to store the cell expression 
        for i in range(M):
            mtx[i]=self.cells[i].E
            print("cluster %s, cell %s"%(self.ID,i))
        mtx=csc_matrix(mtx)  # convert to csc matrix for efficent column operations 
        return mtx
        
    # get the mean, sigma for the expression of cells in the cluster
    def __getAvgVarE(self):
        AE=[]
        pcount=0.01
        R=[] 
        (M,N)=self.EMatrix.shape
        for i in range(N):
            ei=self.EMatrix.getcol(i)
            imu=np.mean(ei)
            ivar=np.mean(ei.power(2))-imu**2
            AE.append(imu)
            R.append(ivar)
            print("cluster %s, gene %s"%(self.ID,i))
        R=[item+pcount for item in R]
        return [AE,R]

# In[5]: Class Path(Edge)
## Class Path (Edge)
class Path:
    def __init__(self,fromNode,toNode,GL,dTD,dTG,fChangeCut=0.6):
        print("build edge %s->%s"%(fromNode.ID,toNode.ID))
        self.fromNode=fromNode                                                  # from Node
        self.toNode=toNode                                                      # to Node
        
        self.FC=self.__getFC(GL)                                                      # fold change for all genes
        self.diffF=self.__getDiffGene(fChangeCut)                                     # get differnetial genes based on log fold change 
        self.diffT=self.__getDiffGeneTest(GL)                                         # get differential genes based on student t- test
        self.diffG=[item for item in self.diffF if item in self.diffT]                # get differnetial genes based on fold change and student t-test
       
        self.ptf=self.fromNode.getExpressedTF(dTD.keys(),GL)                         # expressed TFs
        self.etf=self.__getetf(dTD,dTG,GL,fChangeCut)                            # transcription factors and diff TFs
        
        #---------------------------------------------------------------------- 
        self.B=self.__getTransition(dTD,dTG,GL,fChangeCut)                       # transition offset
        self.Q=self.__getProcessVariance(GL,MU=self.fromNode.E)                      # initial process variance
          
    #--------------------------------------------------------------
    # public functons

    #-----------------------------------------------------------------
    # private functions 
    # calculate fold (log) change between fromNode (cluster) and toNode (cluster)
    def __getFC(self,GL):
        def logfc(x, y):
            return y - x

        AE=self.fromNode.E
        BE=self.toNode.E
        FC=[[abs(logfc(AE[i],BE[i])),logfc(AE[i],BE[i]),i,AE[i],BE[i]] for i in range(len(AE))]

        FC.sort(reverse=True)
        FC=[[GL[item[2]],item[1],item[3],item[4]] for item in FC]
        pdFC=pd.DataFrame(data=FC)
        pdFC.columns=['gene','logfc','A','B']
        pdFC=pdFC.set_index(['gene'])
        return pdFC

    # get differential genes between clusters along the path
    def __getDiffGene(self,FCUT):
        print("get differential genes based on logfc for edge %s->%s"%(self.fromNode.ID,self.toNode.ID))
        DG=[item for item in self.FC.index if abs(self.FC.loc[item]['logfc'])>FCUT]
        return DG

    #-------------------------------------------------------------------
    # get differential genes between clusters along the path
    # using student t-test
    def __getDiffGeneTest(self,GL):
        print("get differential genes based on T test for edge %s->%s"%(self.fromNode.ID,self.toNode.ID))
        cut=5e-2
        
        XM=len(self.fromNode.cells)
        YM=len(self.toNode.cells)
        N=len(GL)
        
        X=lil_matrix((XM,N))
        Y=lil_matrix((YM,N))
        
        for i in range(len(self.fromNode.cells)):
            X[i]=self.fromNode.cells[i].E
            
        for i in range(len(self.toNode.cells)):
            Y[i]=self.toNode.cells[i].E
        
        X=X.tocsc()
        Y=Y.tocsc()
        
        TT=[]
        for i in range(len(GL)):
            Xi=X.getcol(i).toarray()
            Yi=Y.getcol(i).toarray()
            pxy=ttest_ind(Xi,Yi)[-1]
            if pxy<cut:
                TT.append([pxy,GL[i]])
        TT.sort()
        DG=[item[1] for item in TT]
        return DG

    # # get enriched TFs based on significantly diff genes
    #---------------------------------------------------------------
    # dMi: input sequence scanning result
    # n: number of sequences in input
    # dTD: dictionary TF->DNA
    # review erniched TF

    def __getetf(self,dTD,dTG,GL,FCUT):
        # strategy 1 (new): 
        # using manwhiteneyu test => find TF whose target genes are "signiciantly" differential along the edge  (compared with all target genes -background)
        def getEnrichTF():
            pcut=0.1
            K=[item for item in dTD.keys() if item in self.ptf]   # only consider TFs that are expressed in the fromNode (>15% of cells in the node)
            BC=[abs(self.FC.loc[item]['logfc']) for item in GL]
            entf=[]
            for i in K:
                iTargets=[item.upper() for item in dTD[i]]
                iFC=[abs(self.FC.loc[item]['logfc']) for item in iTargets]
                pAB=mannwhitneyu(iFC,BC,alternative="greater")[1]
                if pAB<pcut:
                    entf.append([pAB,i])
            entf=sorted(entf,key=lambda x:x[0])
            return entf
        
        print("infer TFs based on their target genes for edge %s->%s"%(self.fromNode.ID, self.toNode.ID))
        etf=getEnrichTF()
        etf=[item for item in etf if item[1] in self.ptf]
        return etf

    # Lasso regresion model for each path
    def __getTransition(self,dTD,dTG,GL,FCUT=0.6):
        G = self.FC
        etfID = [item[1] for item in self.etf]
        dR={0:2,1:-2,2:0} 
        try:
            [X, Y,U,D] = buildTrain(G, dTG, etfID,GL,FCUT)
            dR = {0: U, 1: D, 2: 0}
            print("filtering TFs using LASSO regression model for edge %s->%s"%(self.fromNode.ID,self.toNode.ID))
            LR = LogisticRegressionCV(penalty='l1', Cs=[1.5, 2, 3, 4, 5], solver="saga", multi_class='auto')
            LR.fit(X, Y)
            CE = LR.coef_
            print("regression...")
            print(CE)
            petf = parseLR(self.etf, CE)
            # ---------------------------------------------------------
            XX = []
            for i in HGL:
                if i in dTG:
                    tfi = dTG[i]
                    xi = [tfi[item] if item in tfi else 0 for item in etfID]
                else:
                    xi = [0] * len(etfID)
                XX.append(xi)
            YY = LR.predict(XX)
            self.etf = petf
        except:
            YY = [0 if G.loc[item]['logfc'] > FCUT else 1 if G.loc[item]['logfc'] < -1 * FCUT else 2 for item in GL]
        YY = [dR[item] for item in YY]
        return YY
    
    # get process noise 
    def __getProcessVariance(self,GL,MU=None):
        # MU : average at time t-1, vector
        # X1: all observation at time point t
        # X2:  all observations at time point t
        N=len(GL)
        Q=[]
        X1=self.fromNode.EMatrix
        X2=self.toNode.EMatrix
        
        if MU==None:
            for i in range(N):
                x1=X1.getcol(i) # all observation for gene i at time t-1
                mui=np.mean(x1)+self.B[i]
                x2=X2.getcol(i)
                v=np.mean(x2.power(2))-mui**2
                Q.append(v)
        else:
            for i in range(N):
                x2=X2.getcol(i)
                mui=MU[i]+self.B[i]
                v=np.mean(x2.power(2))-mui**2
                Q.append(v)
        pcount=0.01
        Q=[item+pcount for item in Q]
        return Q

# In[6]:Graph(Tree)
## class Graph 
class Graph:
    def __init__(self,Cells,tfdna,etfile,GL,Clusters,pagaConnects,rootnodeID=None,fChangeCut=0.6,ncores=None):
        # native graph attributes
        print("initializing graph...")
        self.Cells=Cells
        self.fChangeCut=fChangeCut
        self.etfile=etfile
        self.GL=GL
        self.W=self.getW() # naive  weight for each of the genes 
        
        
        # nodes 
        self.Nodes=self.__buildNodes(Clusters,ncores)
        self.root= self.__guessRoot(pagaConnects,rootnodeID) 
        self.__connectNodes(pagaConnects) # connect the nodes to build the tree
        self.getNodePR()
        self.__adjustRTFs(ncores) # get eTFs for each of hte nodes
        
        # edges
        [self.dTD,self.dTG]=parseTFDNA(tfdna,GL)
        self.Edges=self.__buildEdges(ncores)
        self.Paths = self.__buildPaths()
        
        # likelihood 
        self.llh=None
        
    ##-------------------------------------------------------------------------
    # public functions (can be reached outside of the class to update and process the graph)
        
    # update the graph (tree)
    def updateGraph(self,prRes,ncores=None):
        print("update the graph...")
        
        GL=[item.upper() for item in prRes.var.index]
        prRes.obs['scdiff_cluster']=[item.Label for item in self.Cells]
        prRes.obs=prRes.obs.astype('category')
        
        # update the paga connectivity
        sc.tl.paga(prRes,groups='scdiff_cluster')
        pagaConnects=pd.DataFrame(data=prRes.uns['paga']['connectivities_tree'].toarray())
        
        # update nodes
        self.Nodes=self.__buildNodes(prRes.obs.scdiff_cluster)
        self.root=self.__guessRoot(pagaConnects,self.root.ID)
        self.__connectNodes(pagaConnects)
        self.__adjustRTFs(ncores) 
        self.getNodePR()
        
        # update edges
        self.Edges=self.__buildEdges()
        self.Paths = self.__buildPaths()
        return prRes
        
    # re-assign the cells (assign new cluster labels to all cells)
    def ReAssign(self,ncores=None):
        # re-assign
        print("re-assigning all cells to the tree")
        
        cellParList=[]
        for i in self.Cells:
            for j in self.Nodes:
                pij=(i,j,self.W)
                cellParList.append(pij)
        
        MPRWork=MPR(AssignEachCell,cellParList)
        Res=MPRWork.poolwork()
        del MPRWork
        
        dcellBest={}  # infer the best node for each cell with probability
        for i in Res:
            [icell,iNode,pi]=i
            if icell not in dcellBest:
                dcellBest[icell]=[iNode,pi]
            else:
                if pi>dcellBest[icell][-1]:
                    dcellBest[icell]=[iNode,pi]
        
        newlli=[]
        ract=[]   # record the cells that got re-assigned 
        for i in dcellBest:
            cellID=i
            [NodeID,pi]=dcellBest[cellID]
            cell=[item for item in self.Cells if item.ID==cellID][0]
            Node=[item for item in self.Nodes if item.ID==NodeID][0]
            ract=ract+[cell.ID] if cell.Label!=Node.ID else ract
            # update the cell label 
            cell.Label=Node.ID
            newlli.append(pi)
        print("# of re-assigned cells: %s"%(len(ract)))
        #pdb.set_trace()
        
        # update the node cells
        for i in self.Nodes:
            i.cells=[item for item in self.Cells if item.Label==i.ID]
            
        newlli=sum(newlli)
        return newlli
           
    # get the likelihood for the given assignment  (current node cells)
    def getLikelihood(self,ncores=None):
        print("calculate the likelihood for current Graph cell assignment...")
        llhParList=[]
        for i in self.Cells:
            iNode=[item for item in self.Nodes if item.ID==i.Label][0]
            llhParList.append((i,iNode,self.W))
        
        MPRWork=MPR(getLikelihoodEach,llhParList)
        Tlli=MPRWork.poolwork()
        del MPRWork
        Tlli=sum(Tlli)
        return Tlli
         
    # estimate the prior probability for each cluster (node)
    def getNodePR(self):
        TotalCells=len(self.Cells)
        for i in self.Nodes:
            niCells=len(i.cells)
            i.PR = niCells*1.0/TotalCells
            
    # get W (weight for each of the genes)
    def getW(self,MW=0.5):
        # MW: minimum weight
        W = []
        GL=self.GL
        M=len(self.Cells) # M cells
        N=len(GL)         # N genes
        
        mtx=lil_matrix((M,N)) # lil_matrix (sparse) to store the cell expression 
        for i in range(M):
            mtx[i]=self.Cells[i].E
            print("pull expression from cell %s"%(i))
        
        mtx=csc_matrix(mtx)       #  convert to csr matrix, efficient for column operations
        indptr=mtx.indptr
        print("get gene weights ...")
        for i in range(N):
            inz=indptr[i+1]-indptr[i]
            W.append(inz)
        W=[max(MW,item*1.0/M) for item in W]
        return W
        
    #----------------------------------------------------------------------------------------
    # private functions (can only be reached from other functions in the graph calss)
    # building the path (from root to the leaves) of the tree
    def __buildPaths(self):
        print("building paths...")
        def getCompletePath(en):
            # en: end node
            for i in self.Edges:
                if i.toNode == en:
                    return getCompletePath(i.fromNode) + [i]
            return []

        CP = []  # complete path
        for i in self.Nodes:
            if not i.C:
                cp =getCompletePath(i)
                if cp!=[]:
                    CP.append(cp)
        return CP
        
    # build edges
    def __buildEdges(self,ncores=None):
        print("building edges ...")
        edgeParList=[]
        for i in self.Nodes:
            edgeParList.append((i,self.GL,self.dTD,self.dTG,self.fChangeCut))
            
        MPRWork=MPR(buildEachEdge,edgeParList,ncores)
        P=MPRWork.poolwork()
        del MPRWork  # delete the multi-threading worker to avoid memory leak
        P=[item for item in P if item!=None]
        return P
        
    # get eTF (expression based TFs) for each of the nodes 
    # note, adjustRTF has to stay in Graph class, although it's for each of the nodes
    # the inference requires on the complete graph (e.g., parent-children relationships)
    # the expresson of the TF must be unique (different to the parent, *and* different to at least one sibling node)
    def __adjustRTFs(self,ncores=None):
        print("adjusting RTFs...")
        GL=self.GL
        tflist=self.etfile
        # get RTFs (representating TFs) based on its own expression for each node (the TF expression is different to both parent and siblings)
        tflistpath=pkg_resources.resource_filename(__name__,"tfdata/HumanTFList.txt") if tflist==None else tflist
        try:
            with open(tflistpath,'r') as f:
                TFs=f.readlines()
                TFs=[item.strip().split()[0] for item in TFs]
        except:
            print("error! Please check your input TF List file")
            sys.exit(0)
        
        RTFParList=[]
        eTFs=[item for item in GL if item in TFs]
        for i in self.Nodes:
            RTFParList.append((i,eTFs,GL,self.fChangeCut))
        
        
        MPRWork=MPR(adjustRTFEachNode,RTFParList,ncores)
        Res=MPRWork.poolwork()
        del MPRWork  # delete the multi-threading worker to avoid memory leak
       
        
        for i in Res:
            inode=[item for item in self.Nodes if item.ID==i[1]][0]
            inode.eTFs=i[0]
        
     
    # build nodes
    def __buildNodes(self,Clusters,ncores=None):
        print("building nodes...")
        print("start clustering ...")
        
        for j in self.Cells:
            jID=j.ID
            Yj=Clusters[jID]
            j.Label=int(Yj)
        
        ClusterList=sorted(list(set(Clusters)),key=lambda x:int(x))
        ClusterList=map(lambda x:int(x), ClusterList)
        
        nodeParaList=[]
        for i in ClusterList:
            nodeID=i
            nodeCells=[item for item in self.Cells if item.Label==nodeID]
            nodeParaList.append((nodeID,nodeCells))
        
        MPRWork=MPR(buildEachNode,nodeParaList,ncores)
        AC=MPRWork.poolwork()
        del MPRWork  # delete the multi-threading worker to avoid memory leak
        AC=sorted(AC,key=lambda x:x.T)
        return AC
    
    
    # guess the root node of the tree
    def __guessRoot(self,pagaConnects,rootnodeID=None):
        self.Nodes=sorted(self.Nodes,key=lambda x:x.T)
        if rootnodeID==None:
            root=self.Nodes[0]
        else:
            root=[item for item in self.Nodes if item.ID==int(rootnodeID)][0]
        return root
   
    # connect each node (DFS)
    def __connectS(self,S,Visited,pagaConnects):
        iConnects=pagaConnects.loc[S.ID]+pagaConnects[S.ID]
        Visited.append(S.ID)
        for j in iConnects.index:
            jnode=[item for item in self.Nodes if item.ID==j][0]
            if iConnects[j]>0 and (jnode.ID not in Visited):
                jnode.P=S
                print("%s->%s"%(S.ID,jnode.ID))
                self.__connectS(jnode,Visited,pagaConnects)
                
        
    # connect nondes=> infer the parent node for each of the nodes
    def __connectNodes(self,pagaConnects):
        print("connecting nodes ....")
        
        # delete old parent-child relationship
        for i in self.Nodes:
            i.P=None
            i.C=[]
            
        self.__connectS(self.root,[],pagaConnects)
            
        # add Children node information
        for inode in self.Nodes:
            if inode.P!=None:
                inode.P.C+=[inode]
        

# In[11]: Class MPR
## Multiprocessing running
class MPR:
    def __init__(self,Func,DataList,ncores=None):
        self.Func=Func
        self.DataList=DataList
        self.ncores=ncores
        
    def poolwork(self):
        pool=Pool(processes=self.ncores,maxtasksperchild=1)
        Res=pool.map_async(self.Func,self.DataList)
        pool.close()
        pool.join()
        Res=Res.get()
        return Res

# In[12]: Functions for multi-threading 
## can't use object functions as it will copy the entire object to the memory-> huge mem cost

## build each node
## inodePar is a tuple (nodeID,cells) for each node; cells-> all the cells in the node
def buildEachNode(inodePar):
    (nodeID,nodeCells)=inodePar
    print("building node  %s..."%(nodeID))
    CC = Cluster(nodeCells, int(nodeID))
    return CC

## build each edge
def buildEachEdge(iEdgePar):
    (toNode,GL,dTD,dTG,fChangeCut)=iEdgePar
    if toNode.P:
        p1=Path(toNode.P,toNode,GL,dTD,dTG,fChangeCut)
        return p1
    return None   

##adjustRTFEachNode : infer the eTFs for each of the nodes
def adjustRTFEachNode(RTFParList):
    (Node,eTFs,GL,fcut)=RTFParList
    print("infer expression-based TFs(eTFs) for node %s ..."%(Node.ID))
    if Node.P:
        NodeParent=Node.P
        NodeSib=[item for item in Node.P.C if item!=Node]
        NodeSibCells=[] if NodeSib==[] else [item.cells for item in NodeSib]
        NodeParentCells=NodeParent.cells
        peTFs=[]
        for j in eTFs:
            jdex=GL.index(j)
            [flag,pvp,fcp]=tellDifference(Node.cells,NodeParentCells,NodeSibCells,jdex,fcut)
            if flag:
                peTFs.append([pvp,j,fcp])
        peTFs.sort()
    else:
        peTFs=[]
    return [peTFs,Node.ID]

## assign cell function, private, can only be callced by ReAssign 
## this one has to be public (for multi-threading purpose)
## calculate the assign probability of cell-> Node
def AssignEachCell(cellParList):
    (cell,Node,W)=cellParList
    print("cell : %s"%(cell.ID))
    pi=Node.getAssignProbability(cell,W)
    return [cell.ID,Node.ID,pi]
       

##get the likelihood for node i based on the current cell assignment
#  cell
def getLikelihoodEach(llhParList):
    (cell,Node,W)=llhParList
    print("calculate the likelihood for the cell assignment of cell %s -> node %s"%(cell.ID,Node.ID))
    pi=Node.getAssignProbability(cell,W)
    return pi

# In[7]:Global functions
## 2. Global functions

## tell whether the expression is unique in the specified node
## 1,-1 (unique, 1: higher, -1: lower), 0 (nonunique)
def tellDifference(nodeCells,nodePCells,nodeSibCells,geneIndex,fcut=0.6):
    #print("tell whether current node is unique compared to its parent and siblings for TF %s"%(geneIndex))
    X=[item.E[0,geneIndex] for item in nodeCells]
    XP=[item.E[0,geneIndex] for item in nodePCells]
    fcp=np.mean(X)-np.mean(XP)
    pvp=ranksums(X,XP)[1]
    pcut=0.05
    
    # if no sibling nodes
    if len(nodeSibCells)==0:
        if (pvp<pcut) and (fcp>fcut or fcp<-1*fcut):
            return [1,pvp,fcp]
    
    # if has sibling nodes
    for sNodeCells in nodeSibCells:
        Y=[item.E[0,geneIndex] for item in sNodeCells]
        fcs=np.mean(X)-np.mean(Y)
        pvs=ranksums(X,Y)[1]
        if (pvp<pcut and pvs<pcut) and ((fcp>fcut and fcs>fcut) or (fcp<-1*fcut and fcs<-1*fcut)):
            return [1,pvp,fcp]
        
    return [0,pvp,fcp]

# parse input tfdna file
def parseTFDNA(tfdna,GL):
    RTD=TabFile(tfdna).read('\t') # dictionary to store the TF-DNA info
    DEFAULTACTIVITY=1.0
    try:
        if len(RTD[0])==2:
            TD=[[item[0].upper(),item[1].upper(),DEFAULTACTIVITY] for item in RTD[1:] if len(item)>1]
        elif len(RTD[0])==3:
            TD=[[item[0].upper(),item[1].upper(),float(item[2])] for item in RTD[1:] if len(item)>2]
        else:
            TFs=RTD[0][1:]
            genes=[item[0] for item in RTD[1:]]
            RTDM=[[float(k) for k in item[1:]] for item in RTD[1:]]
            TD=[]
            for i in range(len(genes)):
                for j in range(len(TFs)):
                    TD.append([TFs[j],genes[i],RTDM[i][j]])
    except:
        print("check the format of input TF-DNA interaction file")
        sys.exit(0)

    [dTD,dTG]=getTFDNAInteraction(TD,GL)
    return [dTD,dTG]

# get TF-Gene interactions
def getTFDNAInteraction(TD,GL):
    dTD = {}  # TF->DNA
    dTG = {}  # DNA->TF
    for i in TD:
        if i[2]>0 and i[1].upper() in GL:
            if i[0] not in dTD:
                dTD[i[0]] = [i[1]]
            else:
                
                dTD[i[0]].append(i[1])

            if i[1] not in dTG:
                dTG[i[1]] = {}
                dTG[i[1]][i[0]] = i[2]
            else:
                dTG[i[1]][i[0]] = i[2]
    return [dTD,dTG]

# building traning dataset for regression
def buildTrain(G,dTG,ptf,GL,Fcut=1):
    print("build training set with imbalanced sampling")
    # G: differential genes for a given path
    # dTD: DNA->TF dictionary
    # TF candidate
    Ncut=Fcut/2.0
    UP=[item for item in G if item[1]>Fcut]
    DN=[item for item in G if item[1]<-1*Fcut]
    NN=[item for item in G if abs(item[1])<Ncut]


    U=sum([item[1] for item in UP])/len(UP)
    D=sum([item[1] for item in DN])/len(DN)

    UP=[item[0].upper() for item in UP]
    DN=[item[0].upper() for item in DN]
    NN=[item[0].upper() for item in NN]


    XU=[]
    XD=[]
    XN=[]

    YU=[]
    YD=[]
    YN=[]

    HGL=[item.upper() for item in GL]
    for i in HGL:
        if i in dTG:
            tfi=dTG[i]
            xi=[tfi[item] if item in tfi else 0 for item in ptf]
            if i in UP:
                yi=0
                XU.append(xi)
                YU.append(yi)
            elif i in DN:
                yi=1
                XD.append(xi)
                YD.append(yi)
            elif i in NN:
                yi=2
                XN.append(xi)
                YN.append(yi)

    X=XU+XD+XN
    Y=YU+YD+YN
    
    # to solve the imbalanced training set issue, use over-sampling techqniue- SMOTE
    sm=SMOTE(random_state=0)
    Xs,Ys=sm.fit_sample(X,Y)

    Xs=list(Xs)
    Ys=list(Ys)

    return [Xs,Ys,U,D]

# parse Logistic regression result
def parseLR(etf,LRC):
    LRC=[max(item) for item in LRC.T]
    out_etf=[]
    for i in range(len(LRC)):
        if LRC[i]>0:
            out_etf.append(etf[i]+[LRC[i]])
    return out_etf


# logMessgae
def logMessage(logText,logfile):
    print(logText)
    logfile.write(logText)
    logfile.flush()

# In[8] InferGraph function
## inferGraph based on given inputs
def inferGraph(scg,output,tfdna,tfList,fChangeCut,ncores,rootnodeID,llhcut,MAXLOOP=5):   
    if os.path.exists(output)==False:
        os.mkdir(output)
    logfile=open("%s/runninglog.txt"%(output),'a')
    # log the start time
    tnow="The program starts at : %s \n"%(str(datetime.datetime.now()))
    logMessage(tnow,logfile)
    
    # log the parameters 
    logTextPara="""\
    The program runs with the following arguments
    
    $scdiff2
    -i %s
    -o %s
    -t %s
    --etfListFile %s
    --log2fc %s
    --ncores %s
    --root %s
    --llhCut %s
    --maxloop %s \n\
    \n"""%(scg,output,tfdna,tfList,fChangeCut,ncores,rootnodeID,llhcut,MAXLOOP)
    
    logMessage(logTextPara,logfile)
    
    #Read in the prerun results from the prerun (in h5ad format)
    print("loading back prerun results (h5ad) ...")
    prRes=anndata.read_h5ad(scg)
    
    
    #Convert full matrix to sparse_matrix to reduce the memory usage
    # Expression matrix (sparse)
    prRes.X=csr_matrix(prRes.X)

    # Genes
    GL=[item.upper() for item in prRes.var.index]

    # clusters
    clusters=prRes.obs.leiden

    # paga connectivity
    pagaConnects=pd.DataFrame(data=prRes.uns['paga']['connectivities_tree'].toarray())
    
    # log reading cells
    logText1="reading cells ...\n"
    logMessage(logText1,logfile)
    
    # list to store all cells
    AllCells=[]
    
    for i in range(len(prRes.obs.index)):
        iid=prRes.obs.index[i]
        ti=float(prRes.obs.time[i])
        li=prRes.obs.label[i]
        ei=prRes.X[i,:]
        ci=Cell(iid,ti,ei,li)
        AllCells.append(ci)
        print("load cell: "+str(i))
       
    
    #log clustering
    logText2="clustering cells ...\n"
    logMessage(logText2,logfile)
   
    #load clusters from the prerun results
    clusters=prRes.obs.leiden
    
    
    # log building graph (tree)
    logText3="building graph(tree) ...\n"
    logMessage(logText3,logfile)
    
    G1=Graph(AllCells,tfdna,tfList,GL,clusters,pagaConnects,rootnodeID=rootnodeID,fChangeCut=fChangeCut,ncores=ncores)
   
    #drawing graphs
    if os.path.exists(output)==False:
        os.mkdir(output)

    scg_name=scg.split('/')[-1]
    
    # writing out Graph...
    viz(scg_name,G1,output,prRes)
    
    
    # if MAXLOOP (iterative refinment), calculate the inital llh
    if MAXLOOP>0:
        G1.llh=G1.getLikelihood(ncores) # old likelihood
        ollh=G1.llh
        ILLH=ollh     # initial LLH 
        # log iterative refinment
        logText4="likelihood: %s\n"%(ollh)
        logMessage(logText4,logfile)
   
    
    for loop in range(MAXLOOP):
        logLoopText="->loop: %s \n"%(loop)
        logMessage(logLoopText,logfile)
        
        nllh=G1.ReAssign(ncores)
        G1.llh=nllh
        increase_llh=(nllh-ollh)/abs(ILLH)
        # log iterative refinment
        logText5="likelihood: %s -> likelihood increase this loop: %s\n"%(nllh,increase_llh)
        logMessage(logText5,logfile)
        prRes=G1.updateGraph(prRes,ncores)
        if increase_llh<llhcut:
            break 
        ollh=nllh      # update ollh<-nllh
        
        # update the visualziation file
        
        # log writing visualziation page
        logText6="updating the javascript powered visualization file (%s.html) under the InteractiveViz folder\n"%(scg_name)
        logMessage(logText6,logfile)
        viz(scg_name,G1,output,prRes)
    
    # update the clustering and trajectory plots if there is a PGM iterative refinement
    if 'scdiff_cluster' in prRes.obs:
        logTextPlot="The stopping criteria is met, quit the loop \n\n Updating the PGM refined clustering (UMAP), trajectory (PAGA), and DE genes plots \n"
        logMessage(logTextPlot,logfile)
        sc.settings.figdir = '%s/figures'%(output)
        sc.tl.paga(prRes,groups='scdiff_cluster')
        sc.pl.paga(prRes,show=False,save="_Traj.pdf")
        sc.tl.umap(prRes,init_pos='paga')
        sc.pl.umap(prRes,color=['scdiff_cluster','time'],legend_loc="on data",show=False,save="_clustering.pdf")
        sc.tl.rank_genes_groups(prRes, 'scdiff_cluster', method='wilcoxon')
        sc.pl.rank_genes_groups(prRes, n_genes=25, sharey=False,show=False, save="_global_DE_genes.pdf")
        # update the visualziation with the updated prRes
        viz(scg_name,G1,output,prRes)
    
        # log writing h5ad
        logText7="writing the results to a %s/%s file ...\n"%(output,scg_name)
        logMessage(logText7,logfile)
        prRes.write_h5ad("%s/%s"%(output,scg_name),compression=9)
        
    # log ending 
    logText8="job completed!\n"
    logMessage(logText8,logfile)
    
    tnow="The program ends at : %s \n\n"%(str(datetime.datetime.now()))
    logMessage(tnow,logfile)
    logfile.close()

# In[9]: Main
## main function
def main():
    # parse input arguments
    parser=argparse.ArgumentParser(description="scdiff2 main")
    parser.add_argument('-i','--input',required=True,help='h5ad result from pre-run')
    parser.add_argument('-o','--output',required=True,help='output directory')
    parser.add_argument('-t','--tfdna',required=True, help='TF-DNA interaction data')
    parser.add_argument('--etfListFile',required=False,default=None,help='By default, this program recognizes 1.6k TFs (collected in human and mouse). Users are able ' +
                                                                         'to provide a customized list of TFs　using this option (e.g, for another species).')                                                       
    parser.add_argument('--log2fc',required=False,default=0.6, help='By default, scdiff uses log2 Fold change 0.6(=>2^0.6~=1.5) as the cutoff ' +
                                                                    'for differential genes (together with student t-test p-value cutoff 0.05). ' +
                                                                    'Users can customize this log fold change cutoff.')
    parser.add_argument('--ncores',required=False,default=4, help='# of allocated cpu cores for multi-threading the job (4 by default)')
    parser.add_argument('--root',required=False,default=None, help='Set the root (of the tree) as an input cluster ID (e.g., 0 from the prerun result)')  
    
    parser.add_argument('--llhCut',required=False,default=0.05, help='The convergence likelihood cutoff, the program stops if the cell ' +
                                                                'assignment likelihood improvement is smaller than this cutoff (e.g. 0.05: 5 percent)') 
    parser.add_argument('--maxloop',required=False,default=5, help='The max # of loops allowed for the PGM based iterative refinment. Set it to 0 ' +
                                                                'to directly use the clustering and trajectory results from the prerun program (scanpy based). ' +
                                                                'Only the regulatory networks (TFs) and the interactive visulzation page ' +
                                                                'will be learned and generated')
    
    args = parser.parse_args()
    # input arguments
    scg=args.input
    output=args.output
    tfdna=args.tfdna
    tfList=args.etfListFile
    fChangeCut=float(args.log2fc)
    ncores=int(args.ncores)
    rootNode=args.root
    llhcut=float(args.llhCut)
    MAXLOOP=int(args.maxloop)
    
    inferGraph(scg,output,tfdna,tfList,fChangeCut,ncores,rootNode,llhcut,MAXLOOP)

# In[10]: Program Entry
if __name__=="__main__":
    main()
    

