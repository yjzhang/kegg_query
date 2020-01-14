import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'kegg.db')
SPECIES_MAP = {'human': 'human', 'homo_sapiens': 'human', 'Human': 'human', 'mouse': 'mouse', 'mus_musculus': 'mouse', 'Mouse': 'mouse', 'both': 'all', 'all': 'all'}
BASIC_SPECIES = {'human', 'mouse', 'all'}

def get_all_genes(species='human', db_dir=DB_DIR):
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species != 'all':
        C.execute('SELECT DISTINCT gene FROM genes WHERE species=?', (species,))
    else:
        C.execute('SELECT DISTINCT gene FROM genes')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_all_cells(species='human', db_dir=DB_DIR):
    """
    Returns a list of all unique cell names, or tissue names.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species != 'all':
        C.execute('SELECT DISTINCT pathwayName FROM pathway_genes WHERE species=?', (species,))
    else:
        C.execute('SELECT DISTINCT pathwayName FROM pathway_genes')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_all_cell_cls(species='human', db_dir=DB_DIR):
    """
    Returns a list of (pathwayName, kegg id) tuples.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species != 'all':
        C.execute('SELECT DISTINCT pathwayName, pathwayID FROM pathway_genes WHERE species=?,' (species,))
    else:
        C.execute('SELECT DISTINCT pathwayName, pathwayID FROM pathway_genes')
    results = C.fetchall()
    conn.close()
    return results

def get_all_species(db_dir=DB_DIR):
    return list(BASIC_SPECIES)


@lru_cache(maxsize=None)
def get_cell_genes(cell, species='human', db_dir=DB_DIR):
    """
    Given the name of a cell, this returns a list of all genes associated with that cell.
    Alternatively, if cells_or_tissues is 'tissues', this returns a list of
    all genes associated with that tissue.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species != 'all':
        C.execute('SELECT genes FROM pathway_genes WHERE pathwayName=? AND species=?', (cell, species))
    else:
        C.execute('SELECT genes FROM pathway_genes WHERE pathwayName=?', (cell,))
    results = C.fetchall()
    results = [x[0].split(',') for x in results]
    conn.close()
    return results

def get_all_cells_genes_ids(species='human', db_dir=DB_DIR):
    """
    Returns a list of tuples (pathwayId, pathwayName, [genes]) for all pathways.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species != 'all':
        C.execute('SELECT pathwayID, pathwayName, genes FROM pathway_genes WHERE species=?', (species,))
    else:
        C.execute('SELECT pathwayID, pathwayName, genes FROM pathway_genes')
    results = C.fetchall()
    results = [(x[0], x[1], x[2].split(',')) for x in results]
    conn.close()
    return results


def hypergeometric_test(genes, species='human', return_header=False, db_dir=DB_DIR):
    """
    Uses a hypergeometric test to identify the most relevant cell types.

    Returns:
        list of 4-tuples: id, cell name, p-value, overlapping genes
        in order of ascending p-value.
    """
    if species not in BASIC_SPECIES:
        species = SPECIES_MAP[species]
    from scipy import stats
    genes = [x.upper() for x in genes]
    genes = set(genes)
    all_genes = get_all_genes(species=species)
    all_genes = [x.upper() for x in all_genes]
    all_cells = get_all_cells_genes_ids(species=species, db_dir=db_dir)
    cell_p_vals = []
    # each cell should have 4 items: id, cell type, p-value, overlapping genes
    for cell_id, cell_name, cell_genes in all_cells:
        cell_genes = [x.upper() for x in cell_genes]
        overlapping_genes = list(genes.intersection(cell_genes))
        if len(overlapping_genes) == 0:
            continue
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(cell_genes), len(genes))
        cell_p_vals.append((cell_id, cell_name, 1 - pv, overlapping_genes))
    cell_p_vals.sort(key=lambda x: x[2])
    if return_header:
        header = ['Pathway ID', 'Pathway Name', 'P-value', 'Overlapping Genes']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals
