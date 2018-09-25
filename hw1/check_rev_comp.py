f1=open("lambda.fasta", 'r')
f2=open("lambda.rev.fasta", 'r')

lines1 = f1.readlines()[1:]
lines2 = f2.readlines()[1:]

elements1=0
elements2=0

ellist1 = []
ellist2 = []

for line in lines1:
    new_line = line.strip('\n')
    for element in new_line:
        elements1 +=1
        ellist1.append(element)

for line in lines2:
    new_line = line.strip('\n')
    for element in new_line:
        elements2 +=1
        ellist2.append(element)

print elements1, elements2
#returns 48502 48502
ellist2rev = ellist2[::-1]
print ellist1[245:250]
print ellist2rev[245:250]
#returns complements
