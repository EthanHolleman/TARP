from Transposer.search import Search
from Transposer.process_sams import *

BTI = '/media/ethan/EH_DATA/GMAX_1.1_BTI/GMAX_1.1_BTI'
con_file = '/media/ethan/EH_DATA/TEST_ARGS/Clusters/Gmr3INTACT_cluster_240_consensus'
out_file = '/media/ethan/EH_DATA/TEST_ARGS/Sam_Files/Gmr3INTACT_cluster_240_consensus_intact.sam'
acc = '/media/ethan/EH_DATA/GMax1.1_assembly/chr2acc'
BDB = '/media/ethan/EH_DATA/Gmax1.0_Blast_DB/GMAX_1.0_BDB'
out = '/media/ethan/EH_DATA/TEST_OLD_BTI'
s = Search(BTI=BTI, con_file=con_file, out_file=out_file, acc=acc, BDB=BDB, type="I", num_old_els=500)
s.make_element_set()
print(len(s.element_set))

sort_sam = sort_sams([s])
write_fasta(sort_sam, output=out, name='TEST_SEARCH_METHODS')
