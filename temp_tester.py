from Transposer.search import Search
from Transposer.process_sams import *
'''
BTI = '/media/ethan/EH_DATA/GMAX_1.1_BTI/GMAX_1.1_BTI'
con_file = '/media/ethan/EH_DATA/TEST/Clusters/Gmr3SOLO_cluster_299'
out_file = '/media/ethan/EH_DATA/TEST/Sam_Files/Gmr3SOLO_cluster_299_consensus_intact.sam'
acc = '/media/ethan/EH_DATA/GMax1.1_assembly/chr2acc'
BDB = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/GMAX_1.0_BDB'
out = '/media/ethan/EH_DATA/TEST_OLD_BTI/test.fa'
s = Search(BTI=BTI, con_file=con_file, out_file=out_file, acc=acc, BDB=BDB, type="I", num_old_els=500, intact_len=5000)
print('making set')
s.make_element_set()
print(len(s.element_set))
print('sorting')
sort_sam = sort_sams([s])
print('sorted')
rename = rename_elements(sort_sam)
write_fasta(sort_sam, output=out)


# solo element removal seems to be working well
from Transposer.element import Element
from Transposer.process_sams import prune

a = Element(name='a', accession='b', chr=1, startLocation=100, endLocation=200, length=100, status='S', seq='A')
b = Element(name='a', accession='b', chr=1, startLocation=110, endLocation=210, length=100, status='S', seq='A')
c = Element(name='a', accession='b', chr=1, startLocation=400, endLocation=200, length=100, status='S', seq='A')
d = Element(name='a', accession='b', chr=1, startLocation=500, endLocation=200, length=100, status='S', seq='A')
e = Element(name='a', accession='b', chr=1, startLocation=510, endLocation=200, length=100, status='S', seq='A')

l = [a,b,c,d,e]

g = prune(l)

for i in g:
    print(i.startLocation, 's')



from fasta_tools import read_as_tuples

els = read_as_tuples(fasta_file='/media/ethan/EH_DATA/TARP_Runs/TEST_GMR3_2.1_Clstr_2/Results/TEST_GMR3_2.1_Clstr_2.fa')
hold = []
for h, s in els:
    h = h.split(' ')
    start = int(h[-2][:-1])
    hold.append(Element(h[0], accession=h[0], chr=0, startLocation=int(start), endLocation=0, length=0, status='P', seq=s))

print(len(els), 'element length list')

p = prune(hold, n=100000)

print(len(list(p)), 'prune length')
'''
from Transposer.element import Element
from Clustering.ClstrFile import ClstrFile
from Transposer.process_sams import pruner

a = ClstrFile(path='/media/ethan/EH_DATA/TARP_Runs/TEST_GMR3_2.1_Clstr_2_pruning/Clusters/in_con_clstr.clstr')
print(list(a.get_singles()))

c_path = '/media/ethan/EH_DATA/TARP_Runs/TEST_GMR3_2.1_Clstr_2_oct28/Results/TEST_GMR3_2.1_Clstr_2_oct28.csv'
from Transposer.process_sams import prune
import csv
elements = []
with open(c_path) as c:
    reader = csv.reader(c)
    next(reader)
    for row in reader:
        elements.append(Element(name=row[0], accession=row[1], chr=int(row[2]), startLocation=int(row[3]), endLocation=0, length=0, status='N', seq='A'))
print(len(elements))
a = pruner(elements)

f = open('/media/ethan/EH_DATA/TARP_Runs/TEST_GMR3_2.1_Clstr_2_oct28/Results/TEST_GMR3_2.1_Clstr_2_oct29.csv', 'w')
write = csv.writer(f)
a = ([i.startLocation] for i in a)
for r in a:
    write.writerow(r)
