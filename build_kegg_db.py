import pickle
import sqlite3

with open('pathway_data_human_reduced.pkl', 'rb') as f:
    pathways_human = pickle.load(f)

with open('pathway_data_mouse_reduced.pkl', 'rb') as f:
    pathways_mouse = pickle.load(f)


conn = sqlite3.connect('kegg_query/data/kegg.db')
c = conn.cursor()

c.execute('CREATE TABLE pathway_genes(pathwayID text, pathwayName text, species text, genes text)')
for point in pathways_human:
    if 'genes' not in point or 'id' not in point or 'name' not in point:
        print(point)
        continue
    genes = ','.join(point['genes'])
    c.execute('INSERT INTO pathway_genes VALUES (?, ?, ?, ?)', (point['id'], point['name'], 'human', genes))
for point in pathways_mouse:
    if 'genes' not in point or 'id' not in point or 'name' not in point:
        print(point)
        continue
    genes = ','.join(point['genes'])
    c.execute('INSERT INTO pathway_genes VALUES (?, ?, ?, ?)', (point['id'], point['name'], 'mouse', genes))
c.execute('CREATE INDEX pathway_genes_index ON pathway_genes(pathwayID, species)')

genes_human = set()
for point in pathways_human:
    if 'genes' in point:
        genes_human.update(point['genes'])
genes_mouse = set()
for point in pathways_mouse:
    if 'genes' in point:
        genes_mouse.update(point['genes'])
c.execute('CREATE TABLE genes (gene text, species text)')
for gene in genes_human:
    c.execute('INSERT INTO genes VALUES (?, ?)', (gene, 'human'))
for gene in genes_mouse:
    c.execute('INSERT INTO genes VALUES (?, ?)', (gene, 'mouse'))
c.execute('CREATE INDEX genes_index ON genes(gene)')

conn.commit()
conn.close()
