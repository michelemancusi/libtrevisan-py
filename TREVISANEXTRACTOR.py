"""
Trevisan's extractor
"""


import time
import numpy as np
from random import choices 
import finitefield
from finitefield import *
from bitarray import bitarray


"""
Definizione Funzioni
"""

def GFpow(x,y):                                                                #definizione funzione che fa l'elevazione a potenza come moltiplicazione, perchè, con la libreria finitefield, in GF(2) l'elevazione di 0 a qualunque potenza dava 1          
    if y==0:
        return x**y
    else:
        a=x
        for i in range(y-1):
            x=x*a
        return x            
            
def inttobit(x):                                                               #definizione funzione che trasforma un intero in binario bigendian
    if x>0:
        return [x//2**k % 2 for k in range(int(np.floor(np.log2(x))+1))]
    if x==0:
        return [0]
    
def bittoint(x):                                                               #definizione funzione che converte da binario ad intero
    a=0
    for i in range(len(x)):
        a=a+x[i]*2**i
    return a

def WDcomputeS(i,m,t,t_req):                                                   #definizione funzione weak design
    GFt=GF(t)                                                                  #creazione del campo finito GF(t) (la notazione in GF è bigendian)
    alf=[]
    b=0
    S=[]  
    c=np.ceil(np.log2(m)/np.log2(t_req)-1)  
    mask=(1<<int(np.log2(t)))-1                                                #creazione di una maschera per poter fare operazioni bit a bit
    
    for j in range(int(c)+1):                                                  #preparo i coefficienti del polinomio, devo contare da 0 fino a c compreso.
        alf.insert(0,(i&(mask<<(j*int(np.log2(t)))))>>j*int(np.log2(t)))       #inserisco gli elementi nella lista alf (cioè alpha). Nell'inserire, si procede da destra verso sinistra                                                              
        
    for j in range(int(c)+1):                                                  #trasformo i coefficienti in binario e li faccio diventare elementi di un campo finito GF(t)
        alf[-j-1]=inttobit(alf[-j-1])
        alf[-j-1]=GFt(alf[-j-1])
            
    for a in range(int(t_req)):                                                #calcolo dell'i-esimo set S_i
        b=0
        Sa=0
        a=inttobit(a)                                                          #converto "a" in bit  
        a=GFt(a)                                                               #faccio diventare "a" un elemento di GF(t)               
    
        for j in range(int(c)+1):                                              #calcolo il valore del polinomio p(x) per i vari x, polinomio che qui ho chiamato con "b"      
           b=b+alf[-j-1]*GFpow(a,j)                                       
           
        b=b.coeffs                                                             #crea la lista con gli elementi di "b" tramite .coeffs, partendo da "b" scritto come elemento di GF(t)
        b=bittoint(b)                                                          #trasformo "b" in intero
        Sa=Sa^b                                                                #costruisco per ogni "a" un S^(a)_i, che rappresenta l'elemento dell'i-esimo set S_i determinato dalla coppia (x, p(x)), cioè dalla coppia (a,b). Ora sto inserendo nella metà a destra di Sa il valore in binario di "b"
        a=a.coeffs                                                             #crea la lista con gli elementi di "a" tramite .coeffs, partendo da "a" scritto come elemento di GF(t) 
        a=bittoint(a)                                                          #converto "a" in intero
        Sa=Sa^(a<<int(np.log2(t)))                                             #inserisco "a" scritto in binario in Sa nelle restanti posizioni di Sa, ovvero nella metà di sinistra di Sa    
        S.append(Sa)                                                           #costruisco S inserendo tutti i vari Sa come elementi di S
  
    return S

def OneBitExt(sottoseed, source, error):
    source_bitarray=bitarray(source)                                           #trasformo seed e source in bittarray
    sottoseed_bitarray=bitarray(sottoseed)
    c=[]
    r=0
    n=len(source)
    l=int(np.ceil(np.log2(n)+2*np.log2(2/error)))
    GF2l=GF(2**l)                                                              #creazione del campo finito GF(2^l) (la notazione in GF è bigendian)
    s=int(np.ceil(n/l))
    source=source+[0]*(s*l-n)                                                  #vogliamo dividere la sorgente di lunghezza "n" in "s" sottostringhe tutte di lunghezza "l", quindi allunghiamo la lista source mettendoci (s*l-n) zeri alla fine  
                                                                 
    for i in range(s):                                                         #divido la stringa di input in "s" sottostringhe
        c.append(source[i*l:(i+1)*l])
            
    alf=sottoseed[0:l]                                                         #prendo la prima metà del sottoseed che va dal bit in posizione 0 al bit in posizione l e la inseriso nella variabile alf, poichè la lunghezza "t" del sottoseed è definita come 2*l con l=ceil(np.log2(n)+2*log2(2/error))
    
    for i in range(s):                                                         #Step Reed-Solomon: trasformo le varie sottostringhe in elementi del campo finito GF(2^l)
        c[i]=GF2l(c[i])
       
    alf=GF2l(alf)                                                              #trasformo alf in un elemento del campo finito GF(2^l)
    
    for i in range(1,s+1):                                                     #calcolo il polinomio del Reed-Salomon code: nell'articolo alf è elevato alla (s-i), mentre su internet ho trovato che deve essere elevato alla (i-1); inoltre, sia su internet che nell'articolo, x_i (cioè c[i] in questo caso) parte da x_1, mentre io lo faccio partire da x_0 per come ho costruito le sottostringhe
        r=r+c[i-1]*GFpow(alf,(i-1))   
   
    
    r=r.coeffs                                                                 #crea la lista con gli elementi di "r" tramite .coeffs, partendo da "r" scritto come elemento di GF(2^l)
    r=bitarray(r)
    
    b=0                                                                        #Step Hadamard: calcolo la parità del prodotto bit a bit tra "r" e la seconda metà del sottoseed
    
    for j in range(l):                                                         
        b=b^(sottoseed[j+l]&r[j])
        
    return b        
            
     
"""
Definizione parametri e costanti
"""
start_time = time.time()

r=2*np.e
n=int(input("Qual è la lunghezza della stringa in input?")) 
alpha=float(input("Qual è la min-entropy della stringa in input?")) 
eps=float(input("Qual è l'errore per bit per la costruzione finale?"))

k=alpha*n;
m=int(np.floor((k-4*np.log2(1/eps)-6)/r));

#k=r*m+4*np.log2(1/eps)+6

t_req=int(2*np.ceil(np.log2(n)+2*np.log2(2/eps)))
t=int(2**(np.ceil(np.log2(t_req))))                                            #Trasformiamo t_req nella potenza di 2 più vicina e più grande, poichè il WD funziona con t che sono potenze di 2 quindi sarà t>t_req
d=t**2                         

"""
Inserimento source e seed
"""
q= input("Enter the file name of the SOURCE (with .txt extension): ")
f1=open(q,"r")
a=f1.read()

w= input("Enter the file name of the SEED (with .txt extension): ")
f2=open(w,"r")
b=f2.read()

seed=[]
source=[]

for i in range(d):
    seed.append(int(b[i],2))
    
for i in range(n):
    source.append(int(a[i],2))
    

"""
Algoritmo Trevisan's extractor
"""

rho=[0]*m

for i in range(m):
    print(i/m*100)
    S=WDcomputeS(i,m,t,t_req)
    b=[0]*t_req                          
    
    for j in range(t_req):               
        b[j]=seed[S[j]]
        #print("S[%d]", S[j])
       # print(seed[S[j]])
        
    rho[i]=OneBitExt(b,source,eps)

print()
print("La stringa in uscita casuale è",rho)
print()
print("--- Il programma ha impiegato %s minuti per girare ---" %((time.time() - start_time)/60)) 
print()