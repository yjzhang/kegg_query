# TODO: use the kegg API or some other tool to download the KEGG pathway database for human and mouse.

# use bioservices?
# https://bioservices.readthedocs.io/en/master/kegg_tutorial.html#building-a-histogram-of-all-relations-in-human-pathways

# extract all relations from all pathways
from bioservices.kegg import KEGG
s = KEGG()
s.organism = "hsa"

def simpleParse(data_string):
    """Extract name, genes, id from data_string."""
    point = {}
    mode = None
    for line in data_string.split('\n'):
        if mode == 'gene':
            if not line.startswith(' '):
                mode = None
            else:
                line_data = line.strip().split(';')
                gene = line_data[0].split()[-1]
                point['genes'].append(gene)
        if mode != 'gene':
            if line.startswith('NAME'):
                line0 = line[4:].strip()
                point['name'] = line0
            elif line.startswith('ENTRY'):
                line0 = line[5:].strip().split()[0]
                point['id'] = line0
            elif line.startswith('GENE'):
                mode = 'gene'
                point['genes'] = []
                line_data = line.strip().split('; ')
                gene = line_data[0].split()[-1]
                point['genes'].append(gene)
    return point

# retrieve more than 260 pathways so it takes time
pathway_strings = [s.get(x) for x in s.pathwayIds]
pathway_data = [s.parse(x) for x in pathway_strings]
pathway_points = [simpleParse(x) for x in pathway_strings]

# reducing pathway to only necessary data
#pathway_points = []
#for i, pathway in enumerate(pathway_data):
#   try:
#       point = {'name': pathway['NAME'],
#                'genes': [x.split()[1] for x in pathway['GENE'].keys()],
#                'id': pathway['ENTRY'].split()[0]}
#       pathway_points.append(point)
#   except:
#       print('ERROR: unable to parse pathway string {0}'.format(i))

import pickle
with open('pathway_data_human_reduced.pkl', 'wb') as f:
    pickle.dump(pathway_points, f)
with open('pathway_data_human_raw.pkl', 'wb') as f:
    pickle.dump(pathway_strings, f)



s = KEGG()
s.organism = "mmu"

# retrieve more than 260 pathways so it takes time
pathway_strings = [s.get(x) for x in s.pathwayIds]
pathway_data = [s.parse(x) for x in pathway_strings]

# reducing pathway to only necessary data
pathway_points = [simpleParse(x) for x in pathway_strings]

import pickle
with open('pathway_data_mouse_reduced.pkl', 'wb') as f:
    pickle.dump(pathway_points, f)
with open('pathway_data_mouse_raw.pkl', 'wb') as f:
    pickle.dump(pathway_strings, f)
