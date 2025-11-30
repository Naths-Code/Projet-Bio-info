print ("Hello World ! The project is about to begin ! ")

import os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Toujours configurer l'email pour NCBI
Entrez.email = "votre.email@exemple.com"

def charger_donnees_utilisateur():
    """
    Fonction interactive pour récupérer des séquences.
    Retourne une LISTE d'objets SeqRecord prêts à être traités.
    """
    print("\n--- MENU D'IMPORTATION ---")
    print("1. Télécharger via Numéro d'Accession (NCBI)")
    print("2. Charger depuis un fichier FASTA local")
    print("3. Coller une séquence brute")
    
    choix = input("Votre choix (1/2/3) : ").strip()

    sequences_recuperees = []

    # --- CAS 1 : ACCESSION ---
    if choix == "1":
        acc_id = input("Entrez le numéro d'accession (ex: M62320) : ").strip()
        print(f"Connexion à NCBI pour {acc_id}...")
        try:
            # On demande le format GenBank pour avoir les annotations si besoin, ou fasta simple
            handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text")
            sequences_recuperees = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            print("Téléchargement réussi !")
        except Exception as e:
            print(f"Erreur lors du téléchargement : {e}")

    # --- CAS 2 : FICHIER LOCAL ---
    elif choix == "2":
        chemin = input("Entrez le chemin du fichier FASTA (ex: data/mon_fichier.fasta) : ").strip()
        if os.path.exists(chemin):
            try:
                sequences_recuperees = list(SeqIO.parse(chemin, "fasta"))
                print(f"{len(sequences_recuperees)} séquences chargées.")
            except Exception as e:
                print(f"Erreur de lecture du fichier : {e}")
        else:
            print("Erreur : Le fichier n'existe pas.")

    # --- CAS 3 : SÉQUENCE BRUTE ---
    elif choix == "3":
        seq_str = input("Collez votre séquence (ATGC...) : ").strip().upper()
        nom = input("Donnez un nom à cette séquence : ").strip()
        # Création manuelle d'un objet SeqRecord (comme le ferait Biopython)
        mon_record = SeqRecord(
            Seq(seq_str),
            id=nom,
            description="Séquence importée manuellement"
        )
        sequences_recuperees.append(mon_record)

    else:
        print("Choix invalide.")

    return sequences_recuperees

# --- TEST DU MODULE ---
if __name__ == "__main__":
    # Ce bloc ne s'exécute que si on lance ce script directement
    seqs = charger_donnees_utilisateur()
    
    if seqs:
        print("\nRésumé des données chargées :")
        for s in seqs:
            print(f"- ID: {s.id} | Longueur: {len(s.seq)}")
            
            # ICI, on peut proposer à l'utilisateur de couper
            # C'est là qu'on réutilise votre logique validée :
            # troncon = s.seq[debut-1 : fin]