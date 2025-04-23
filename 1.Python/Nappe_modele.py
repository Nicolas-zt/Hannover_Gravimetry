import flopy
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from scipy.ndimage import zoom

# Grille simple
nlay = 2
nrow = 600
ncol = 600
delr = delc = 10.0

with rasterio.open("../4.Data/MNT_TOULON.vrt") as src:
    mnt = src.read(1)

   
def MNT():
    
    # Charger une carte d'altitude par raster
     
        
        
    plt.figure(figsize=(8, 6))
    plt.imshow(mnt, cmap="terrain",vmin = -50, vmax = 300)
    plt.colorbar(label="Altitude (m)")
    plt.title("Modèle Numérique de Terrain (MNT)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
    


# Adapter à ta grille MODFLOW






def modflow():
    
    
    # top = np.interp(np.linspace(0, mnt.shape[0], nrow),
    #                 np.arange(mnt.shape[0]),
    #                 np.mean(mnt, axis=1))[:, None] * np.ones((nrow, ncol))
    
    # Calcul des facteurs de réduction
    zoom_factors = (nrow / mnt.shape[0], ncol / mnt.shape[1])
    
    # Redimensionnement bilinéaire
    top = zoom(mnt, zoom_factors, order=1)  # order=1 = interpolation bilinéaire

    botm1 = top - 50  # sol d'épaisseur 50m
    botm2 = top -100 # nappe d'épaisseur 50m
    
    botm = [botm1,botm2]
    # Création du modèle
    modelname = "nappe_perte"
    mf = flopy.modflow.Modflow(modelname, exe_name="mf2005")
    

    
    dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc, top=top, botm=botm)
    
    # Conditions de base
    ibound = np.ones((nlay, nrow, ncol), dtype=int)
    strt = 5.0 * np.ones((nlay, nrow, ncol))
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
    
    # Paramètres hydrogéologiques
    hk = 10.0
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk)
    
    # Recharge générale
    rch = flopy.modflow.ModflowRch(mf, rech=0.001)
    
    # 💡 Définition d'une perte d'eau locale (puits fictif)
    # Exemple : case au centre (50, 50) de la couche 0
    # wel_spd = {0: [[0, 50, 50, -500]]}  # débit négatif = extraction (m3/jour)
    
    #  wel = flopy.modflow.ModflowWel(mf, stress_period_data=wel_spd)
    
    # Solveur et options de sortie
    pcg = flopy.modflow.ModflowPcg(mf)
    oc = flopy.modflow.ModflowOc(
        mf, stress_period_data={(0, 0): ['print head', 'save head', 'save budget']}, compact=True
    )
    
    # Exécution
    mf.write_input()
    success, buff = mf.run_model(silent=True)
    
    # Lecture des résultats
    headobj = flopy.utils.HeadFile(f"{modelname}.hds")
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])
    
    # Visualisation
    plt.figure(figsize=(8, 6))
    plt.imshow(head[0], cmap="Blues", extent=(0, ncol * delr, 0, nrow * delc))
    plt.colorbar(label="Niveau piézométrique (m)")
    plt.title("Effet d'une perte d'eau locale (puits fictif)")
    plt.xlabel("Distance Est (m)")
    plt.ylabel("Distance Nord (m)")
    plt.grid(False)
    plt.show()
    
    return head





if __name__ == "__main__":
    
   modflow()