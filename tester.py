from Clustering.ClstrFile import ClstrFile

a = ClstrFile('/media/ethan/EH_DATA/TARP_Runs/Big_Run/GMR_111/Clusters/Gmr111SOLO.clstr')
for clstr in a.clusters_set:
    print(clstr.num_elements)
a.write_cluster_fastas(original_fasta_path='/media/ethan/EH_DATA/TARP_Runs/Big_Run/GMR_111/Clusters/Gmr111SOLO', path='/media/ethan/EH_DATA/TARP_Runs/Big_Run/GMR_111/Clusters')

for clstr in a.clusters_set:
    print(clstr.fasta)

print(a.path)
