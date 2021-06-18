import numpy as np



class Mutable_seq:
    def __init__(self,seq,nam):
        self.seq=seq
        self.length=len(seq)
        self.name=nam
        #self.history
        
        
#    def get_len(self):
 #       return len(self)
    
    
    def mutate(self,prob):
        ct=0
    #    length = len(mutable_seq)
        NucList = {'A':['C','G','T'],'C':['A','G','T'],'G':['A','C','T'],'T':['A','C','G']}
        vals=np.random.binomial(size=self.length,n=1,p=prob)
        for i, nt in enumerate(self):
            if vals[i]:
                nt=nt.upper()
                self.seq = "".join( [ "".join(self.split()[:i]) , random.choice(NucList[nt]) , "".join(self.split()[i+1:])])
                ct=ct+1
        
        print(self.length, ct)
        return self.upper()

    
    
    def invert(self):
        return self[::-1].upper()
    
    
    
    
    
    
    
    #s at start
    #l #r, left third right third
    #m middle
    #e at end


    def insert(self, mutated_seq, pos):
        length=self.length
        convert={'s':0,'e': length-1, 'l': round(length/3), 'r':round(2*length/3),'m':round(length/2)  }


        if type(pos)==str:
            poss = convert[pos]
        elif type(pos)==int:
            poss=pos
        self = "".join( [ "".join(self[:poss]) , mutated_seq , "".join(self[poss:])])
        return self.upper()
    
    
    
    
    def fastify(seq, chrom,start,end):
        return "\n".join(     [ "/".join([chrom,str(start),str(end)])   , seq    ]   )