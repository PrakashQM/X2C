#! /usr/bin/python
import re
import copy
import sys

if len(sys.argv) == 1:
    print "To_uncontract.py.  Use:  To_uncontract.py <contracted_basis_file> [<uncontracted_basis_file>]"
    exit()

input_basis = ""
if len(sys.argv) == 2:
    input_basis = sys.argv[1]
    split_input_basis = input_basis.split('.')
    output_basis = ".".join(split_input_basis[:-1]) + "-un." +  split_input_basis[-1]

def find_searchquery():
    searchquery=[]
    f1 = open(input_basis,'r')
    lines = f1.readlines()
    for  i in range(0, len(lines)):
            if re.search('^spherical',lines[i]):
                print 'sperical'
                print ' '
            if re.search('^cartesian',lines[i]):
                print 'cartesian'
                print ' '    
            if re.search('^\!',lines[i]):
                print lines[i],
            if re.search('^[A-Z].?     0',lines[i]):
                searchquery.append(lines[i])
    return searchquery
    
def print_uncontracted_basis(searchquery):
    k=0
    basis_holder={}
    with open(input_basis) as f1:
            lines = f1.readlines()
            for  i in range(0, len(lines)):
                line=lines[i]          
                if line.startswith(searchquery):
                    for  j in range(i, len(lines)):
                        pv=lines[j]
                        if (pv[0]=="*") :
                            break
                        else:
                            basis_holder[k]=pv.split()[0:]
                            k= k+1
                                   
    seen=copy.copy(basis_holder)
    ll=0
    for i in range(len(basis_holder)):
        if (i==0):
            print '****'
            print basis_holder[0][0],' ' ,basis_holder[0][1]
        else:
            if re.search('^(S|P|D|F|G|H|I)',basis_holder[i][0]):                   
                    if  not basis_holder[i]==seen[ll]:
                        pv = basis_holder[i][1]
                        Ipv=int(pv)            
                        for k in range(1,Ipv+1):
                                print basis_holder[i][0],' ', 1, ' ',1.00
                                print ' ', basis_holder[i+k][0], ' ',1.00
                    ll=i 
                                       
sys.stdout=open(output_basis,"w+")

print "! Uncontracting the basis set"
searchquery=list()
searchquery=find_searchquery()
print " "
print " "
print " "
for j in range(len(searchquery)):
    print_uncontracted_basis(searchquery[j])                 
print '****'
 
 

                    
                
        

        


