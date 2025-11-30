import sys
import warnings
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# === CONFIGURATION GLOBALE ===
# On filtre les warnings pour ne pas polluer la console
warnings.filterwarnings("ignore")

# Email obligatoire pour s'identifier aupres du NCBI
Entrez.email = "nathandissezpro@gmail.com"

# Liste officielle demandee par la consigne du projet
LISTE_OFFICIELLE = [
    "M62320",   # VIH-1 Groupe M Sous-type A
    "K03455",   # VIH-1 Groupe M Sous-type B
    "U46016",   # VIH-1 Groupe M Sous-type C
    "K03454",   # VIH-1 Groupe M Sous-type D
    "AJ006022", # VIH-1 Groupe N
    "L20587",   # VIH-1 Groupe O
    "M15390",   # VIH-2 (ROD)
    "M33262",   # SIVmac (Macaque)
    "X52154"    # SIVcpz (Chimpanze)
]

def extraire_genes_interessants(record_genbank):
    """
    Analyse une fiche GenBank pour trouver et extraire les sequences ENV et TAT.
    Retourne un dictionnaire contenant les sequences proteiques traduites.
    """
    genes_trouves = {"env_prot": None, "tat_prot": None}
    
    # On parcourt chaque annotation du fichier GenBank
    for feature in record_genbank.features:
        # On ne s'interesse qu'aux regions codantes (CDS)
        if feature.type == "CDS":
            
            # Recuperation securisee du nom du gene (parfois dans 'gene', parfois dans 'product')
            # On met tout en minuscule pour faciliter la comparaison
            nom_gene = feature.qualifiers.get("gene", [""])[0].lower()
            produit = feature.qualifiers.get("product", [""])[0].lower()
            
            # --- Detection de l'Enveloppe (Env) ---
            if "env" in nom_gene or "envelope" in produit:
                # On verifie si le NCBI fournit deja la traduction en proteine
                if "translation" in feature.qualifiers:
                    genes_trouves["env_prot"] = feature.qualifiers["translation"][0]
            
            # --- Detection de Tat ---
            # Note : Tat est souvent eparpille (epissage).
            # On fait confiance a l'annotation 'translation' du NCBI qui recolle les morceaux sans introns.
            elif "tat" in nom_gene:
                if "translation" in feature.qualifiers:
                    genes_trouves["tat_prot"] = feature.qualifiers["translation"][0]
    return genes_trouves
1
def pipeline_principal():
    print("\n=== PIPELINE D'ACQUISITION ET D'EXTRACTION DE DONNEES ===")
    
    # === ETAPE 1 : CHOIX DES DONNEES CIBLES ===
    print("\nQuelles sequences voulez-vous utiliser ?")
    print("1. Utiliser la liste officielle du sujet (VIH-1, VIH-2, SIV...)")
    print("2. Saisir manuellement des numeros d'accession")
    
    choix = input("Input 1 ou 2) : ").strip()
    cibles = []

    if choix == "1":
        print(f"[INFO] Chargement de la liste officielle ({len(LISTE_OFFICIELLE)} virus).")
        cibles = LISTE_OFFICIELLE
    elif choix == "2":
        entree = input("Entrez les ID separes par des virgules (ex: M62320, K03455) : ")
        # On nettoie la saisie en enlevant les espaces inutiles)
        cibles = [x.strip() for x in entree.split(",") if x.strip()]
    else:
        print("[ERREUR] Choix non reconnu. Arret du programme.")
        return

    # Containers pour stocker nos resultats
    sequences_env = []
    sequences_tat = []
    
    print(f"\n[INFO] Demarrage du traitement pour {len(cibles)} souches...\n")

    # === ETAPE 2 : BOUCLE PRINCIPALE TELECHARGEMENT et EXTRACTION ===
    for acc_id in cibles:
        try:
            # A. Telechargement via API NCBI
            # On demande le format 'gb' (GenBank) et non 'fasta', necessaire pour avoir les coordonnees des genes (features).
            print(f"> Traitement de {acc_id}...")
            handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # B. Extraction intelligente (Appel de notre fonction dediee)
            resultats = extraire_genes_interessants(record)
            
            # C. Stockage des resultats positifs
            nom_virus = f"{acc_id}" # L'ID devient le nom de la sequence
            
            # Gestion Env
            if resultats["env_prot"]:
                # On cree un objet SeqRecord propre pour Biopython
                seq_obj = SeqRecord(Seq(resultats["env_prot"]), id=nom_virus, description="Env_Protein")
                sequences_env.append(seq_obj)
                print(f"  [SUCCES] Gene Env recupere.")
                print (sequences_env)
            else:
                print(f"  [ATTENTION] Pas de gene Env trouve.")

            # Gestion Tat
            if resultats["tat_prot"]:
                seq_obj = SeqRecord(Seq(resultats["tat_prot"]), id=nom_virus, description="Tat_Protein")
                sequences_tat.append(seq_obj)
                print(f"  [SUCCES] Gene Tat recupere.")
                print (sequences_tat)
            else:
                print(f"  [ATTENTION] Pas de gene Tat trouve.")
            
        except Exception as e:
            print(f"  [ERREUR] Impossible de traiter {acc_id} : {e}")

    # === ETAPE 3 : SAUVEGARDE DES FICHIERS POUR L'ALIGNEMENT ===
    print("\n=== GENERATION DES FICHIERS DE SORTIE ===")
    
    # On n'ecrit le fichier que si on a trouve des sequences
    if sequences_env:
        nom_fichier_env = "env_all.fasta"
        SeqIO.write(sequences_env, nom_fichier_env, "fasta")
        print(f"[FICHIER CREE] '{nom_fichier_env}' contient {len(sequences_env)} sequences.")
    else:
        print("[AVERTISSEMENT] Aucune sequence Env n'a ete trouvee.")
    
    if sequences_tat:
        nom_fichier_tat = "tat_all.fasta"
        SeqIO.write(sequences_tat, nom_fichier_tat, "fasta")
        print(f"[FICHIER CREE] '{nom_fichier_tat}' contient {len(sequences_tat)} sequences.")
    else:
        print("[AVERTISSEMENT] Aucune sequence Tat n'a ete trouvee.")

    print("\n[INFO] Pipeline termine. Prets pour l'alignement Clustal Omega.")

if __name__ == "__main__":
    pipeline_principal()