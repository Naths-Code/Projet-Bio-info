import os
from Bio import Entrez, SeqIO

# 1. Configuration (Mets ton email, c'est obligatoire pour NCBI)
Entrez.email = "ton.email@exemple.com" 

# Dictionnaire des souches cibles
TARGETS = {
    "HIV1_A_U455": "M62320",
    "HIV1_B_HXB2": "K03455",
    "HIV1_C_ETH2220": "U46016",
    "HIV1_D_ELI": "K03454",
    "HIV1_N_YBF30": "AJ006022",
    "HIV1_O_ANT70": "L20587",
    "HIV2_ROD": "M15390",
    "SIVmac_239": "M33262",
    "SIVcpz_GAB1": "X52154"
}

def fetch_and_extract():
    if not os.path.exists("data"):
        os.makedirs("data")

    print(f"Téléchargement de {len(TARGETS)} génomes depuis NCBI...")
    
    # Révision : Utilisation des listes et dictionnaires
    ids = list(TARGETS.values())
    
    # Appel à l'API NCBI (efetch)
    # on demande des fichiers type "gb" (GenBank) qui contiennent les annotations (positions des gènes)
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")

    sequences_found = {"env_prot": [], "tat_prot": [], "env_dna": []}

    for record in records:
        # Retrouver le nom lisible à partir de l'ID (ex: M62320 -> HIV1_A_U455)
        name = [k for k, v in TARGETS.items() if v == record.name or v in record.id][0]
        print(f"Traitement de {name}...")

        # Révision : Parsing d'objets complexes (les Features GenBank)
        for feature in record.features:
            if feature.type == "CDS": # Coding Sequence
                # On cherche le qualifieur "gene" ou "product" qui contient "env" ou "tat"
                gene_name = feature.qualifiers.get("gene", [""])[0].lower()
                product = feature.qualifiers.get("product", [""])[0].lower()

                # Extraction ENV
                if "env" in gene_name or "envelope" in product:
                    # On évite les doublons (parfois env est coupé en gp120/gp41 dans les annotations)
                    # Pour faire simple ici, on prend le premier grand fragment trouvé
                    if len(feature.location) > 2000: 
                        sequences_found["env_dna"].append(f">{name}\n{feature.extract(record.seq)}")
                        # Biopython traduit automatiquement si la traduction est dispo, sinon on traduit
                        prot = feature.qualifiers.get("translation", [feature.extract(record.seq).translate()])[0]
                        sequences_found["env_prot"].append(f">{name}\n{prot}")
                
                # Extraction TAT
                elif "tat" in gene_name:
                    # Tat est souvent épissé (spliced) en 2 exons. Biopython gère ça avec feature.extract()
                    prot = feature.qualifiers.get("translation", [feature.extract(record.seq).translate()])[0]
                    sequences_found["tat_prot"].append(f">{name}\n{prot}")

    # Sauvegarde dans des fichiers FASTA
    for key, fasta_list in sequences_found.items():
        filename = f"data/{key}.fasta"
        with open(filename, "w") as f:
            f.write("\n".join(fasta_list))
        print(f"Fichier créé : {filename} ({len(fasta_list)} séquences)")

if __name__ == "__main__":
    fetch_and_extract()