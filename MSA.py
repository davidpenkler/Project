#MSA.py
# David Penkler
import Functions
option = ""
file_name = ""
Functions.menu(file_name)
while option.upper() != "Q" and option != "7":
    option = raw_input().upper()
    if option == "":
            Functions.menu(file_name)    
    elif option == "O" or option == "1":
        data, file_name, count = Functions.open_MSA()
    elif option == "A" or option == "2":
        if file_name != "":
            Functions.MSA_info(data,file_name,count)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)"
    elif option == "E" or option == "3":
        if file_name != "":
            Functions.MSA_exp(file_name,count,data)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)"  
    elif option == "I" or option == "4":
        if file_name != "":
            Functions.seq_info(data)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)" 
    elif option == "S" or option == "5":
        if file_name != "":
            Functions.Gly_Print(data)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)"
    elif option == "X" or option == "6":
        if file_name != "":
            Functions.fasta_exp(data)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)"
    elif option == "V" or option == "8":
        if file_name != "":
            Functions.MSA_gap_score(data)
        else:
            print "ERROR: An MSA file needs to be opened first (Option 1)"            
    elif option != "Q" and option != "7":
        print "Please enter a valid option (1-7)"
print "Thank you for using Multiple Alignment Analyzer"
        