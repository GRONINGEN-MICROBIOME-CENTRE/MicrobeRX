'''
Hello, dear user!

We are glad you are interested in our awesome project. 
However, we regret to inform you that we are still working on this part of the code. 
Please don't be mad at us, we are doing our best to deliver a high-quality product as soon as possible. 
In the meantime, you can check out our other features or watch some funny cat videos to pass the time. 
We appreciate your patience and understanding. Thank you for your support and stay tuned for updates!

Sincerely,
Angel Ruiz
'''


from io import StringIO
from multiprocessing import Pool
import pandas as pd

from tqdm.notebook import tqdm
import requests
from io import StringIO


ORGANISMS_PATH='data/STRAINS_uniqueSpecies.tsv'

class GetUniprotData:
    
    def __init__(self, query:str, query_type:str):
        
        self.organisms_database = pd.read_csv(ORGANISMS_PATH,sep='\t')
        
        def __generate_chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield lst[i:i + n]
                
        species_chunks=list(__generate_chunks(self.organisms_database.NCBI_taxon_ID.unique(),200))

        data = [f'({query_type}:{query})+AND+taxonomy_id:({"+OR+".join(map(str,d))})' for d in species_chunks]

        self.fields={'accession':'ACCESION',
                'id':'ENTRY_id',
                'reviewed':'STATUS', 
                'annotation_score':'SCORE',
                'protein_existence':'EXISTENCE',
                'protein_name':'PROTEIN_name',
                'cc_function':'FUNCTION',
                'gene_primary':'GENE_name',
                'xref_geneid':'GENE_id',
                'xref_refseq' : 'REFSEQ',
                'organism_name':'ORGANISM_name',
                'organism_id':'ORGANISM_id',
                'lineage' : 'LINEAGE',
                'ec':'EC',
                'go_p':'GO_process',
                'go_c':'GO_component',
                'go_f':'GO_function',
                'protein_families':'PROTEIN_families',
                'ft_domain':'FEATURE_domain',
                'ft_motif':'FEATURE_motif',
                'cc_domain':'COMMENT_domain',
                'ft_topo_dom' : 'DOMAIN_topological',
                'xref_pdb':'PDB',
                'xref_kegg':'KEGG',
                'xref_biocyc':'BIOCYC',
                'xref_interpro':'INTERPRO',
                'xref_pfam':'PFAM',
                'ft_act_site':'FEATURE_active_site',
                'ft_site':'FEATURE_site',
                'ft_binding':'FEATURE_binding',
                'cc_catalytic_activity':'CATALYTIC_activity',
                'rhea':'RHEA',
                'length':'LENGTH',
                'sequence':'SEQUENCE'}

        self.urls=[f"https://rest.uniprot.org/uniprotkb/search?fields={','.join(self.fields.keys())}&format=tsv&query={u}&size=500" for u in data]
        
    def run_search(self):
        results=[]
        for url in tqdm(self.urls):
            r=requests.get(url)
            results.append(r.text)
            '''
            if 'next' in r.links.keys():
                repeat=True
            else:
                repeat=False


            while repeat:
                r=requests.get(r.links['next']['url'])
                results.append(r.text)
                if 'next' in r.links.keys():
                    repeat=True
                else:
                    repeat=False
            '''
        
        self.uniprot_data = pd.concat([pd.read_csv(StringIO(result), sep='\t') for result in results])
        self.uniprot_data.drop_duplicates(inplace=True,ignore_index=True)
        self.uniprot_data.columns=self.fields.values()
        self.uniprot_data.sort_values(by='SCORE',ascending=False,ignore_index=True,inplace=True)