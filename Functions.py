#Functions.py
# David Penkler

# Displays methods in menu format
def menu(file_name):
    print "*"*70
    print "* MULTIPLE ALIGNMENT ANALYZER","*".rjust(40)
    print "*"*70
    print "* Select an option from below:", "*".rjust(39)
    print "*","*".rjust(68)
    print "*\t1) Open a Multiple Alignment File\t(O)","*".rjust(18)
    print "*\t2) Alignment Information\t\t(A)","*".rjust(18)
    print "*\t3) Alignment Explorer\t\t\t(E)","*".rjust(18)
    print "*\t4) Information per Sequence\t\t(I)","*".rjust(18)
    print "*\t5) Analysis of Glycosylation Signatures\t(S)","*".rjust(18)
    print "*\t6) Export to Fasta\t\t\t(X)","*".rjust(18)
    print "*\t7) Exit\t\t\t\t\t(Q)","*".rjust(18)
    print "*","*".rjust(68)
    print "*\t8) View MSA and calculate gap penalties\t(V)","*".rjust(18)
    print "*","*".rjust(68)
    print "*","File:".rjust(50),file_name,"*".rjust(16-len(file_name))
    print "*"*70
    return"Enter option: "

# 1) Loads the contents of MSA.aln file into memory
import re
def open_MSA():
    path = raw_input("Enter a valid path for a Multiple Sequence Alignment: ")
    if ".aln" in path:
        try:
            f = open(path,"r")
            head = f.readline()
            data = f.readlines()
        
            #this nex line addition is crucial for regular expression m to work every time as it terminates at a new line. some clustal files ie clustalw have a nex line at the end of the file
            data[len(data)-1] = data[len(data)-1] + "\n"  
            flag = True
            count = 0
            if head[:7] != "CLUSTAL":
                print "WARNING! This file is not in CLUSTALw format"
            if not data[0].isspace():
                print "WARNING! This file is not in CLUSTALw format"
            
            i = 0
            while data[i].isspace():
                i = i + 1
         
            c = "([^ ]+)\s+([ATGC-]+)\n"
            cobj = re.compile(c)
            m = "([^ATGC-]+)\n"
            mobj = re.compile(m)
            seq_id1 = cobj.findall(data[i])
            count = 1
            pos = i
            seq_lst = []
            seq_lst.append([seq_id1[0][0],""])
            length = len(data[i])-(len(seq_id1[0][1])+1)
            i=i+1
            
            
            while i < len(data):
                if  data[i].count("*") == 0:
                    seq_id = cobj.findall(data[i])
                    if not data[i].isspace() and seq_id[0][0] != seq_lst[0][0]:
                        seq_lst.append([seq_id[0][0],""])            
                        count = count + 1  
                    else:
                        break    
                i = i +1
                
        
            seq_lst.append([" "*len(seq_lst[0][0]),""])
            
            while pos < len(data)-(count):
                for k in range(0,count+1):
                    if k != count:
                        seq_id = cobj.findall(data[pos+k])
                        seq_lst[k][1] = seq_lst[k][1] + seq_id[0][1]
                    else:
                        seq_id = mobj.findall(data[pos+(count)])
                        seq_lst[k][1] = seq_lst[k][1] + seq_id[0][length:]
                pos = pos + (count + 2)
            print "The file",path,"has been successfully loaded"
                  
            return seq_lst, path, count
        except IOError:
            print "ERROR: Input file not found"
            return "","",""
    else:
        print "ERROR: Invalid input file. File must be .aln"
        return "","",""



# 2) Displays MSA information
def MSA_info(data, file_name,count):
    print "\nFile name:\t", file_name,"\n"
    
    seqs = ""
    i = 0
    while i < count:
        seqs = seqs + "\t\t" + data[i][0] + "\n"
        i = i + 1
    print "Sequences:\n", seqs
    
    print "Length: \t", len(data[0][1]), "bp\n"
    
    i = 0
    tlen = 0
    while i < count:
        for k in data[i][1]:
            if k != "-":
                tlen = tlen + 1
        i = i + 1
    print "Average length:\t","%.2f"%(float(tlen)/count),"\n"
    
    match100 = 0
    for k in data[count][1]:
        if k == "*":
            match100 = match100 + 1
    print "% of Matches:\t", "%.2f"%((float(match100)/len(data[0][1]))*100),"\n"
    
    num_match = []
    for k in range(0,count):
        num_match.append(0)
    
    for k in range (0,len(data[0][1])):
        if data[count][1][k] == " ":
            temp = []
            for i in range(0,count-1):
                temp.append(data[i][1][k])
            i = 0
            while i < len(temp):
                num = -1
                for j in range(0,len(temp)):
                    if (temp[i] != "" or temp[i] != "-") and temp[i] == temp[j]:
                        num = num + 1
                        if num != 0:
                            temp[j] = ""
                if num > 0:
                    temp[i] = ""
                    num = num + 1
                    num_match[num] = num_match[num] + 1          
                i = i + 1
    num = 0            
    for k in data[count][1]:
        if k == "*":
            num = num +1
    print "Number of X matches:\n"
    for k in range(2,count):
        print "\t\t[",k,"] = ",num_match[k],"\n"
    print "\t\t[",count,"] = ", num,"\n"
    print "Press [enter] to display the menu again: "                
    return ""

# 3) Alignment segment explorer
def MSA_exp(file_name,count,data):
    start = input("Enter the start of segment: ")
    end = input("Enter the end of segment: ")
    ranges = "[" + str(start) + "-" + str(end) + "]"
    if end-start > 0 and end <= len(data[0][1]):
        print "\nCLUSTAL segment",ranges, "of the", file_name,"alignment\n"
    
        while end - start >= 60:
            for k in data:
                print k[0],"\t\t",k[1][start-1:start+59]
            start = start+60
        if end - start > 0:
            for k in data:
                print k[0],"\t\t",k[1][start-1:end]  
        print "Press [enter] to display the menu again: "  
    else:
        print "ERROR: segment range invalid or out of range"
    return ""

# 4) Displays Sequence information by Sequence ID
def seq_info(data):
    seq_id = raw_input("Enter the sequence ID: ")
    flag = False
    for k in data:
        if k[0] == seq_id:
            flag = True
            print "\nID:\t", seq_id
            length = 0
            ungapped = ""
            bases = [0,0,0,0]
            for i in k[1]:
                if i != "-":
                    ungapped = ungapped + i
                    length = length + 1
                    if i == "A":
                        bases[0] = bases[0] + 1
                    elif i == "T":
                        bases[1] = bases[1] + 1
                    elif i == "G":
                        bases[2] = bases[2] + 1
                    else: 
                        bases[3] = bases[3] + 1
            print "\nLength:", length
            print "\nFrequences per base:"
            print "\n\t\t    [A]:",bases[0],"\n\t\t    [T]:",bases[1],"\n\t\t    [G]:",bases[2],"\n\t\t    [C]:",bases[3]
            print "\nSequence:\n"
            
            start = 0
            while length - start > 60:
                print ungapped[start:start+60]
                start = start+60
            if length - start > 0:
                print ungapped[start:]
    
    if flag == False:
        print "Invalid sequence ID entered"
    print "\nPress [enter] to display the menu again: " 
    return

# Translates a DNA sequence to Amino Acid sequence
def translate(seq):
    #Amino Acid Codon Dictinary
    codon_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"|","TAG":"|","TGA":"|","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G","TGG":"W"}
    
    
    protein_seq = ""
    k = 0
        
    while k <= len(seq)-3:
        codon = ""
        codon = seq[k]+seq[k+1]+seq[k+2]
        protein_seq = protein_seq + codon_dict[codon]
        k = k + 3
    return protein_seq

# Removes the gaps for a proein or DNA sequence
def ungap(data):
    protein_seqs = []
    
    for k in range(0,len(data)-1):
        ungapped = ""
        for i in data[k][1]:
            if i != "-" :
                ungapped = ungapped + i
        protein_seqs.append(ungapped)
                
    return protein_seqs   

# Finds Glycosylation signitures in a protein sequence    
def Gly_Sign(data):
    protein_lst = []
    protein_seqs = ungap(data) 
    for k in range(0,len(data)-1):
        protein_seqs[k] = translate(protein_seqs[k]).lower()
        protein_lst.append([])
    
    for k in range(0,len(protein_seqs)):
        for i in protein_seqs[k]:
            protein_lst[k].append(i)    
    
    protein_seqs = []    
    for k in range(0,len(protein_lst)):
        protein_seqs.append([])
        sequence = ""
        for i in range(0,len(protein_lst[k])):
            if  i < len(protein_lst[k])-2 and protein_lst[k][i] == "k":
                if protein_lst[k][i+1] != "p"and protein_lst[k][i+1] != "|":
                    if protein_lst[k][i+2] == "s" or protein_lst[k][i+2] == "t":
                        protein_lst[k][i] = protein_lst[k][i].upper()
                        protein_lst[k][i+1] = protein_lst[k][i+1].upper()
                        protein_lst[k][i+2] = protein_lst[k][i+2].upper() 
            sequence = sequence + protein_lst[k][i]
        protein_seqs[k] = sequence
    return protein_seqs

# 5) Shows the Glycosylation signitures in a DNA sequence
def Gly_Print(data): 
    protein_seqs = Gly_Sign(data)
    pos = []
    neg = []
    
    for k in range(0,len(protein_seqs)):
        if protein_seqs[k].islower() == False:
            string = protein_seqs[k] + "\n"
            pos.append(string)
        else:
            neg.append(data[k][0])
        
    if len(pos) > 0:
        print "\nGlycosylation signitures found in:\n" 
        for k in range(0,len(pos)):
            start = 0
            print data[k][0]
            while len(pos[k])-start > 60:
                print pos[k][start:start+60]
                start = start +60
            if len(pos[k])-start > 0:
                print pos[k][start:]
        
    if len(neg) > 0 and len(neg) <= 2:    
        print "There are no signatures in:\n",
        print neg[0],"and",neg[1]
    elif len(neg) > 0:
        print "There are no signatures in:\n",
        for k in range(0,len(neg)-2):
            print neg[k],",", 
        print neg[k+1],
        print "and",neg[len(neg)-1]
            
    print "\nPress [enter] to display the menu again: "
    return protein_seqs

# 6) Exports sequences in MSA to either one file or multiple files, with or with out glycosylation signitures
def fasta_exp(data):
    ftype = raw_input("\nIf you want to create multiple files(one per sequence) press M for a single file containing all the sequences press S: ")
    if  ftype.upper() != "S" and ftype.upper() != "M":
        print "Error Please select one of the two option provided"
        print "\nPress [enter] to display the menu again: "
    else:
        base = raw_input("\nDo you want to save the DNA sequence (input D) or the protein sequence (input P): ")    
        if base.upper() != "P" and base.upper() != "D":
            print "Error Please select one of the two option provided"
            print "\nPress [enter] to display the menu again: "
        else: 
            if ftype.upper() == "S":
                if base.upper() == "P":
                    Glydata = raw_input("\nDo you want to mark by case the bases involved in Glycosylation signatures (Y/N): ")
                    path = raw_input("\nPlease input the name of the file to save without extention (.fasta default): ")+".fasta"
                    f = open(path,"w")
                    protein = Gly_Sign(data)
                    if Glydata.upper() == "Y":
                        start = 0
                        for k in range(0,len(protein)):
                            f.write(">"+data[k][0]+"\n")
                            while len(protein[k]) - start > 60:
                                f.write(protein[k][start:start+60]+"\n")
                                start = start + 60
                            if len(protein[k]) - start > 0:
                                f.write(protein[k][start:]+"\n\n")
                        print "\nFile",path, "successfully saved."
                    elif Glydata.upper() == "N":                
                        protein_lst = ungap(data)
                        for k in range(0,len(protein_lst)):
                            start = 0
                            f.write(">"+data[k][0]+"\n")
                            string = translate(protein_lst[k])
                            while len(protein_lst[k])-start > 60:
                                f.write(string[start:start+60]+"\n")
                                start = start + 60
                            if len(protein_lst[k])-start > 0:
                                f.write(string[start:]+"\n\n")
                        print "\nFile",path, "successfully saved."
                elif base.upper() == "D":  
                    dna_lst = ungap(data)
                    path = raw_input("\nPlease input the name of the file to save without extention (.fasta default): ")+".fasta"
                    f = open(path,"w")            
                    for k in range(0,len(dna_lst)):
                        start = 0
                        f.write(">"+data[k][0]+"\n")
                        string = dna_lst[k]
                        while len(dna_lst[k])-start > 60:
                            f.write(string[start:start+60]+"\n")
                            start = start + 60
                        if len(dna_lst[k])-start > 0:
                            f.write(string[start:]+"\n\n")
                    print "\nFile",path, "successfully saved."
            elif ftype.upper() == "M":
                if base.upper() == "P":
                    Glydata = raw_input("\nDo you want to mark by case the bases involved in Glycosylation signatures (Y/N): ")
                    if Glydata.upper() == "Y":
                        start = 0
                        protein = Gly_Sign(data)
                        for k in range(0,len(protein)):
                            f = open(data[k][0].translate(None,",/?\:|")+"_Protein_Glycos.fasta","w")
                            f.write(">"+data[k][0]+"\n")
                            while len(protein[k]) - start > 60:
                                f.write(protein[k][start:start+60]+"\n")
                                start = start + 60
                            if len(protein[k]) - start > 0:
                                f.write(protein[k][start:]+"\n\n")
                            f.close()
                        print "\nProtein files were successfully saved" 
                    elif Glydata.upper() == "N":
                        protein_lst = ungap(data)
                        for k in range(0,len(protein_lst)):
                            start = 0
                            f = open(data[k][0].translate(None,",/?\:|")+"_Protein.fasta","w")
                            f.write(">"+data[k][0]+"\n")
                            string = translate(protein_lst[k])
                            while len(protein_lst[k])-start > 60:
                                f.write(string[start:start+60]+"\n")
                                start = start + 60
                            if len(protein_lst[k])-start > 0:
                                f.write(string[start:]+"\n\n")
                        print "\nPotein files successfully saved."
                elif base.upper() == "D":
                    dna_lst = ungap(data)
                    for k in range(0,len(dna_lst)):
                        start = 0
                        f = open(data[k][0].translate(None,",/?\:|")+"_DNA.fasta","w")
                        f.write(">"+data[k][0]+"\n")
                        string = dna_lst[k]
                        while len(dna_lst[k])-start > 60:
                            f.write(string[start:start+60]+"\n")
                            start = start + 60
                        if len(dna_lst[k])-start > 0:
                            f.write(string[start:]+"\n\n")
                    print "\nDNA files successfully saved."                      
            else:
                print "Error: Please follow the required inputs correctly"
            print "\nPress [enter] to display the menu again: "
    return ""


# 8) A bonus feture that displays the MSA and calculates the total gap penalty for each sequence in the alignment using and opeing score of -10 and an extension score of -1
def MSA_gap_score(data):
    scores = []
    for k in range(0,len(data)-1):
        i = 0
        gap = 0
        while i < len(data[k][1])-1:
            if data[k][1][i] == "-":
                gap = gap + 10
                while data[k][1][i+1] == "-":
                    gap = gap + 1
                    i = i +1
            i = i + 1
                
        scores.append(gap)
    
    print "\nMultiple Sequence alignment:\n"
    start = 0
    while len(data[0][1])-start > 60:
        for k in data:
            print k[0],"\t\t",k[1][start:start+60]
        start = start + 60
    if len(data[0][1])-start > 0:
        for k in data:
            print k[0],"\t\t",k[1][start:]
            
    print "\nSequence Gap Penalties (opening = -10 ; extention = -1): "
    for k in range(0,len(scores)):
        string = "\t\t" + data[k][0] + ": -" + str(scores[k])
        print string
    print "\nPress [enter] to display the menu again: "
        
                   
    return ""

