from kegg_query import hypergeometric_test

genes = """Atp5e
Rps24
Ybx1
Selenow
Rps14
Rps19
Nsa2
Rplp1
Rpl36
Fam81a
Bag1
Gm10073
Rps3a1
Rps16
Rpl35a
Rpl32
Rpl26
Rpl37a
Pfdn5
Rpl17"""
genes = genes.split()

print(hypergeometric_test(genes))
print(hypergeometric_test(genes, species='mouse'))
print(hypergeometric_test(genes, species='all'))
