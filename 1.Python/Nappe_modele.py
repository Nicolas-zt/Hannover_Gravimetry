import flopy
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from scipy.ndimage import zoom, gaussian_filter
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
# from flopy.mf6.modflow.mfgwfrch import ModflowGwfRch
# from flopy.utils.binaryfile import MF6HeadFile

# Grid
nlay = 8 # Number of sub-soil layers
nrow = 600 # Number of rows
ncol = 600 # Nuber of columns
delr = delc = 10.0 # Size of rows / columns in meters


# Time periods
nper = 8 # Number of time periods
perlen = [30,60,45,45,60,30,60,35] # Duration of each period in days
nstp = [3]*8 # Number ofsteps per period, determines the temporal precision



# Opening the dtm file as .vrt
with rasterio.open("../4.Data/MNT/MNT_TOULON.vrt") as src:
    dtm = src.read(1)   

   
def DTM():
    '''
    Creates a 2D plot of a Digital Terrain Model

    '''
    
     
        
        
    plt.figure(figsize=(8, 6))
    plt.contourf(dtm,levels = 15, cmap="terrain")
    plt.colorbar(label="Altitude (m)")
    plt.title("Modèle Numérique de Terrain (MNT)")
    plt.gca().invert_yaxis()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
    

 
# Adapter à ta grille MODFLOW

def layer_creator(pos,size,shape,thickness,depth = 50):
    '''
    Creates a surface to be used with modflow as a bottom layer. the grid is manually adapted 
    to the one used with the modflow model.
    
    Parameters:
        
        - pos: Center(s) of the layer.
        - size: size of the aquifere.
        - shape: shape of the aquifere (ellipse, gaussian, rectangle or flat).
                 
        - thickness: thickness of the aquifere whithin the layer. Different from the the layer's depth
        
    Out:
        
        - layer: matrix of depths at each location of a grid
    '''
    
    layer = np.zeros((nrow, ncol))
    
    x = np.arange(0.5, ncol) * delr
    y = np.arange(0.5, nrow) * delc
    
    X, Y = np.meshgrid(x, y)
    
    if shape == "ellipse":
        
        cx, cy = pos[0], pos[1]  # centers
        rx, ry = size[0], size[1]  # radius
        
        layer = -thickness*np.exp(-(((X - cx)**2 / rx**2 + (Y - cy)**2 / ry**2))) - depth
        
        
    elif shape == "gaussian":
        
        cx, cy = pos[0], pos[1] 
        sx, sy = size[0], size[1]
        
        layer = thickness * np.exp(-(((X - cx)**2)/(2 * sx**2) + ((Y - cy)**2)/(2 * sy**2))) - depth
        
    elif shape == "rectangle":
            
        
        cx, cy = pos[0], pos[1]
        sx, sy = size[0], size[1]
        
        mask = (np.abs(X - cx) <= sx / 2) & (np.abs(Y - cy) <= sy / 2)
        layer[mask] = thickness
        
        layer -= depth
        
    elif shape == "flat":
        
        layer -= depth
    
        pass
    
    
    return layer + np.linspace(0, 5, nrow)[:, None] 

def layer_plotter(layer):
    '''
    Creates a 2D plot of a layer's surface with matplotlib

    '''
    x = np.arange(0, ncol * delr, delr)  # 0 à 6000 m
    y = np.arange(0, nrow * delc, delc)  # 0 à 6000 m
    X, Y = np.meshgrid(x, y)
    
    plt.figure(figsize=(8,6))
    plt.contourf(X,Y,layer,levels = 100, cmap="terrain")
    plt.colorbar(label="Altitude du fond de nappe (m)")
    plt.gca().invert_yaxis()
    plt.title("Surface basse de la nappe (modélisée)")
    plt.show()

def modflow():
    '''
    Runs a region model with different created sub-soils and a real Digital Terrain Model

    '''
    # top = np.interp(np.linspace(0, mnt.shape[0], nrow),
    #                 np.arange(mnt.shape[0]),
    #                 np.mean(mnt, axis=1))[:, None] * np.ones((nrow, ncol))
    
    # Calcul des facteurs de réduction
    zoom_factors = (nrow / dtm.shape[0], ncol / dtm.shape[1])
    
    # Redimensionnement bilinéaire
    top = zoom(dtm, zoom_factors, order=1)  # order=1 = interpolation bilinéaire
    
    
    

    ######################################### Layers ###################################################

    ### Ground thickness
    botm1 = top + layer_creator([0,0],[0,0],'flat',0,50)  # sol d'épaisseur 50m
    
     
    ### 1st layer - water
    
    
    botm2 = botm1 + layer_creator([2000,1000],[1800,1000],"ellipse",500)
    
    ### 2nd layer - water 
    
    botm3 = botm2 + layer_creator([1000,4000],[3000,1500],'ellipse',100)
    
    ### 3rd layer - soil  
    
    botm4 = botm3 + layer_creator([0,0],[0,0],'flat',0,200)
    
    ### 4th layer - water
    
    botm5 = botm4 + layer_creator([3000,3000],[1000,750],'gaussian',-100)
    
    ### 5th layer - soil 
    
    botm6 = botm5 + layer_creator([0,0],[0,0],'flat',0,300)
    
    ### 6th layer - water
    
    botm7 = botm6 + layer_creator([4000,3000],[3000,3000],'ellipse',500)
    
    ###7th layer - soil 
    
    botm8 = botm7 + layer_creator([0,0],[0,0],'flat',0,1000)
    
    ######################################### Layers ###################################################
   
    
    
    
    
    
    ######################################### Modflow setup ############################################
    
    

    
    botm = [botm1,botm2,botm3,botm4,botm5,botm6,botm7,botm8] # array containing each layer
    
    
    ### Création du modèle
    modelname = "nappe_perte"
    mf = flopy.mf6.MFSimulation(sim_name=modelname, version="mf6", exe_name="mf6", sim_ws="./modflow6_model")

    

    tsmult = [1.0] * nper
    
    tdis = flopy.mf6.ModflowTdis(
    mf,
    nper=nper,
    time_units="DAYS",
    perioddata=[(perlen[i], nstp[i], tsmult[i]) for i in range(nper)],
    )
    

    
    gwf = flopy.mf6.ModflowGwf(mf, modelname=modelname)
    
    dis = flopy.mf6.ModflowGwfdis(
    gwf,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
    )
    
    
    
    
    ### Start conditions
    
    #Sea
    # sea_mask = dtm <= 0
    ibound = np.ones((nlay, nrow, ncol), dtype=int)
    # for k in range(nlay):
    #     ibound[k][sea_mask] = 0 
     
    # chd_data = []
    # for row in range(nrow):
    #     for col in range(ncol):
    #         if sea_mask[row, col]:
    #             chd_data.append([0, row, col, 0.0, 0.0])  # couche 0, hleft = hright = 0
    # chd = flopy.modflow.ModflowChd(mf, stress_period_data={0: chd_data})
    
    # Start levels
    strt = 5.0 * np.ones((nlay, nrow, ncol))
    
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    
    
    ### Hydrogeological parameters
    hk_values = [ 
       1E-5, # Ground
       1E-2, # 1st layer
       5E-2, # 2nd layer
       1E-3, # 3rd layer
       1E-6, # 4th layer soil
       1E-3, # 5th layer
       5E-5, # 6th layer soil
       5E-4, # 7th layer
       1E-10 # 8th layer soil
        ]
        
    vk_values = [val * 0.1 for val in hk_values]
    sy_values = [
        0.15, # Ground
        0.25, # 1st layer
        0.30, # 2nd layer
        0.15, # 3rd layer
        0.02, # 4th layer soil
        0.20, # 5th layer
        0.01, # 6th layer soil
        0.10, # 7th layer
        0.005 # 8th layer soil
        ]
    
    hk = np.empty((nlay,nrow,ncol))
    vk = np.empty((nlay,nrow,ncol))
    sy = np.empty((nlay,nrow,ncol))
    
    for i in range(nlay):
        hk[i, :, :] = hk_values[i]
        vk[i, :, :] = vk_values[i]
        sy[i, :, :] = vk_values[i]
    
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk, k33=vk)
    
    sto = flopy.mf6.ModflowGwfsto(
    gwf,
    sy=sy,
    ss=1e-5,
    iconvert=1,
    steady_state={0: False},
    transient={0: True}
    )
    
    ### general reload
    # rch = flopy.modflow.ModflowRch(mf, rech=0.001)
    
    ### Water variations
    
    # Data
    
    
    perioddata = {
        0 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        1 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        2 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        3 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        4 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        5 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0) 
            ],
        6 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        7 : [(0, "RATE", -50.0),
             (1, "RATE", -50.0),
             (2, "RATE", -50.0)
            ],
        }
    
    
    
    
    maw_pkg = flopy.mf6.ModflowGwfmaw(
    gwf,
    pname="MAW",
    boundnames=True,
    print_input=True,
    print_flows=True,
    save_flows=True,
    mover=False,
    packagedata=[
        (0, 0.1, np.min(botm7), 90.0, "SKIN", 7),  # (maw_id, radius, bottom, type, pump elev, head)
        (1, 0.1, np.min(botm3), 90.0, "SKIN", 2),
        (2, 0.1, np.min(botm5), 90.0, "SKIN", 4)
    ],
    connectiondata=[
        (0, 0, (0, 500, 200), np.min(top)  , np.min(botm1), hk_values[0] - 1E-7, 0.2),  # 1st layer
        (0, 1, (1, 500, 200), np.min(botm1), np.min(botm2), hk_values[1] - 1E-7, 0.2),  # 2nd layer
        (0, 2, (2, 500, 200), np.min(botm2), np.min(botm3), hk_values[2] - 1E-7, 0.2),  # 3rd layer
        (0, 3, (3, 500, 200), np.min(botm3), np.min(botm4), hk_values[3] - 1E-7, 0.2),  # 4th layer
        (0, 4, (4, 500, 200), np.min(botm4), np.min(botm5), hk_values[4] - 1E-7, 0.2),  # 5th layer
        (0, 5, (5, 500, 200), np.min(botm5), np.min(botm6), hk_values[5] - 1E-7, 0.2),  # 6th layer
        (0, 6, (6, 500, 200), np.min(botm6), np.min(botm7), hk_values[6] - 1E-7, 0.2),  # 7th layer
        (1, 0, (0, 100, 400), np.min(top)  , np.min(botm1), hk_values[0] - 1E-7, 0.2),
        (1, 1, (1, 100, 400), np.min(botm1), np.min(botm2), hk_values[1] - 1E-7, 0.2),
        (2, 0, (0, 300, 300), np.min(top)  , np.min(botm1), hk_values[0] - 1E-7, 0.2),
        (2, 1, (1, 300, 300), np.min(botm1), np.min(botm2), hk_values[1] - 1E-7, 0.2),
        (2, 2, (2, 300, 300), np.min(botm2), np.min(botm3), hk_values[2] - 1E-7, 0.2),
        (2, 3, (3, 300, 300), np.min(botm3), np.min(botm4), hk_values[3] - 1E-7, 0.2)    
    ],
    perioddata=perioddata
    )
    
   
    
   
    
   
    
   
    # Reload
    
    rch = flopy.mf6.modflow.ModflowGwfrch(
    gwf,
    stress_period_data={
                        0: [ botm1,float(0.001)],
                        1: [ botm1,float(0.001)],
                        2: [ botm1,float(0.001)],
                        3: [ botm1,float(0.001)],
                        4: [ botm1,float(0.001)],
                        5: [ botm1,float(0.001)],
                        6: [ botm1,float(0.001)],
                        7: [ botm1,float(0.001)]
    },  
    pname="RCH-1",
    filename="model.rch"
    )

    
    
    ### Solveur et options de sortie
    
    ims = flopy.mf6.ModflowIms(mf, print_option="SUMMARY", complexity="SIMPLE")
    mf.register_ims_package(ims, [gwf.name])
    
    oc = flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord=f"{modelname}.hds",
    budget_filerecord=f"{modelname}.cbc",
    saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]
    )
    
    
    ######################################### Modflow setup ############################################
    
    
    
    ### Exécution
    mf.write_simulation()
    success, buff = mf.run_simulation()
    
    ### Lecture des résultats
    # headobj = MF6HeadFile("nappe_perte.gwf.nappe_perte.hds")  # ou adapte le chemin si tu as changé le nom du modèle
    headobj = flopy.utils.HeadFile('./modflow6_model/nappe_perte.hds') 
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])
    
    ### Piezometric levels visualisation
    # plt.figure(figsize=(8, 6))
    # plt.imshow(head[0], cmap="Blues", extent=(0, ncol * delr, 0, nrow * delc), origin ='lower')
    # plt.colorbar(label="Niveau piézométrique (m)")
    # plt.title("Effet d'une perte d'eau locale (puits fictif)")
    # plt.xlabel("Distance Est (m)")
    # plt.ylabel("Distance Nord (m)")
    # plt.grid(False)
    # plt.show()
    
    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')
    
    # X, Y = np.meshgrid(np.arange(ncol), np.arange(nrow))
    
    # # Exemple : top et botm[1] (première couche)
    # ax.plot_surface(X, Y, top, cmap='terrain', alpha=0.7)
    # ax.plot_surface(X, Y, botm[0], cmap='Greens', alpha=0.5)
    # ax.plot_surface(X, Y, botm[1], cmap='Blues', alpha=0.4)
    
    # ax.set_title("Représentation 3D des couches du modèle")
    # ax.set_xlabel("Colonne")
    # ax.set_ylabel("Ligne")
    # ax.set_zlabel("Altitude (m)")
    # plt.show()
    
    
    # Ouvrir le fichier de têtes
    try:
        
        headfile = flopy.utils.HeadFile('./modflow6_model/nappe_perte.hds')  # adapter chemin si besoin
        times = headfile.get_times()  # Liste des temps totaux
        
        # Récupérer la tête moyenne à chaque temps
        avg_heads = []
        for t in times:
            head = headfile.get_data(totim=t)
            avg_heads.append(np.nanmean(head))  # Moyenne en ignorant les NaN
        
        # Tracer
        plt.figure(figsize=(8, 4))
        plt.plot(times, avg_heads, marker='o')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Niveau d'eau moyen (m)")
        plt.title("Évolution du niveau d'eau moyen")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
    except :
        pass

    ################################### Area 3D visualisation ##########################################
    
    # Define the grid's size to be the same as the modflow's one
    x = np.arange(0, ncol * delr, delr)
    y = np.arange(0, nrow * delc, delc)
    X, Y = np.meshgrid(x, y)
    
    plotter = pv.Plotter()
    colors = ['lightgreen', 'blue', 'cyan', 'orange', 'blue','orange','blue','orange']
    
    for i in range(len(botm)):
        # Determine the top and bottom surfaces
        z_top = top if i == 0 else botm[i - 1]
        z_botm = botm[i]
    
        # Creates the 3D points for top and bottom surface of each layer
        z1 = z_top.ravel()
        z2 = z_botm.ravel()
        points = np.concatenate([
            np.c_[Y.ravel(), X.ravel(), z1],
            np.c_[Y.ravel(), X.ravel(), z2]
        ])
    
        # 3D grid
        grid = pv.StructuredGrid()
        grid.points = points
        grid.dimensions = (ncol, nrow, 2)
    
        # Adds volume to the plot
        plotter.add_mesh(grid, color=colors[i % len(colors)], opacity=1, show_edges=False)
    
    plotter.add_axes()
    plotter.show()
    

    ################################### Area 3D visualisation ##########################################
    
    # return head

def surface_nappes(choix):
    
    # if choix == "cos" :
        
    #     x = np.linspace(0, 1, ncol)
    #     y = np.linspace(0, 1, nrow)
    #     X, Y = np.meshgrid(x, y)
        
    #     # Surface basse : variation douce (ex : fond de nappe)
    #     botm_surface = -20 + 2 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
    
    
    
    # if choix == "irregulier":
        
    #     # Bruit aléatoire
    #     noise = np.random.normal(0, 1, size=(nrow, ncol))
        
    #     # Lissage pour rendre la surface réaliste
    #     botm_surface = gaussian_filter(noise, sigma=10)
        
    #     # Redimensionner pour l'échelle désirée (par exemple entre -50m et -30m)
    #     botm_surface = (botm_surface - botm_surface.min()) / (botm_surface.max() - botm_surface.min())
    #     botm_surface = -50 + 20 * botm_surface
        
    # else :
        
    #     # Pente douce
    #     pente = np.linspace(0, -20, nrow)[:, None] * np.ones((nrow, ncol))
        
    #     # Bruit aléatoire
    #     noise = np.random.normal(0, 1, size=(nrow, ncol))
    #     noise = gaussian_filter(noise, sigma=8)
        
    #     # Surface finale
    #     botm_surface = pente + 2 * noise
    
    fond_global = -50 + np.linspace(0, 5, nrow)[:, None]  # pente nord-sud
        
    x = np.arange(ncol)
    y = np.arange(nrow)
    X, Y = np.meshgrid(x, y)
    
    # Centre et taille de la nappe locale
    cx, cy = 250, 450  # centre
    rx, ry = 200, 100  # rayons
    
    # Masque ellipsoïdal
    mask = ((X - cx)**2 / rx**2 + (Y - cy)**2 / ry**2) <= 1
    
    
    # Surface locale aléatoire (ondulée)
    local_surface = -40 + gaussian_filter(np.random.normal(0, 1, size=(nrow, ncol)), sigma=10)
    
    # Fusion : on remplace uniquement à l'intérieur du masque
    fond_final = np.where(mask, local_surface, fond_global)
    
    # Optionnel : lisser un peu toute la surface pour transition douce
    fond_final = gaussian_filter(fond_final, sigma=2)
        
    plt.figure(figsize=(8,6))
    plt.contourf(fond_final, cmap="terrain")
    plt.colorbar(label="Altitude du fond de nappe (m)")
    plt.title("Surface basse de la nappe (modélisée)")
    plt.show()




if __name__ == "__main__":
    
   # layer_plotter(layer_creator([2000,1000],[1800,1000],"ellipse",5))
   # layer_plotter(layer_creator([1000,4000],[3000,1500],'ellipse',5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   # layer_plotter(layer_creator([3000,3000],[1000,750],'gaussian',-5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   # layer_plotter(layer_creator([4000,3000],[3000,3000],'ellipse',5))
   # layer_plotter(layer_creator([0,0],[0,0],'flat',0,5))
   modflow()              