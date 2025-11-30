import os
import subprocess

def lancer_clustal_portable(nom_entree, nom_sortie):
    """
    Fonction : Pilote Clustal Omega en allant le chercher dans son dossier specifique.
    """
    # 1. ANCRAGE (Point de depart : dossier du script)
    dossier_script = os.path.dirname(os.path.abspath(__file__))
    
    # 2. CHEMIN VERS L'EXECUTABLE (La modification est ici)
    # On remonte de 'Programme' vers la racine (..), puis on descend dans le dossier Clustal
    dossier_clustal = os.path.abspath(os.path.join(dossier_script, "..", "clustal-omega-1.2.2-win64"))
    executable = os.path.join(dossier_clustal, "clustalo.exe")
    
    # 3. CHEMIN VERS LES DONNEES
    dossier_donnees = os.path.abspath(os.path.join(dossier_script, "..", "Donnees"))
    chemin_in = os.path.join(dossier_donnees, nom_entree)
    chemin_out = os.path.join(dossier_donnees, nom_sortie)

    # 4. VERIFICATIONS
    if not os.path.exists(executable):
        print(f"[ERREUR] Executable introuvable.")
        print(f"Cherché ici : {executable}")
        print("Verifiez que le nom du dossier 'clustal-omega-1.2.2-win64' est exact.")
        return

    if not os.path.exists(chemin_in):
        print(f"[ERREUR] Fichier de donnees introuvable : {chemin_in}")
        return

    print(f"[INFO] Lancement de l'alignement pour : {nom_entree}")

    # 5. EXECUTION
    # On utilise -v pour avoir des details en cas d'erreur
    commande = [
        executable,
        "-i", chemin_in,
        "-o", chemin_out,
        "--force",
        "--outfmt=fasta",
        "-v"
    ]

    try:
        # On capture la sortie
        resultat = subprocess.run(commande, capture_output=True, text=True)

        if resultat.returncode == 0:
            print(f"[SUCCES] Fichier genere : {nom_sortie}")
        else:
            print(f"[ERREUR] Clustal a echoue (Code {resultat.returncode}).")
            # Si stderr est vide, l'erreur est souvent dans stdout
            msg_erreur = resultat.stderr if resultat.stderr else resultat.stdout
            print(f"Detail : {msg_erreur}")

    except Exception as e:
        print(f"[ERREUR TECHNIQUE] {e}")

if __name__ == "__main__":
    lancer_clustal_portable("env_all.fasta", "env_aligne.fasta")
    lancer_clustal_portable("tat_all.fasta", "tat_aligne.fasta")