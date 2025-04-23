import flopy
import numpy as np
import matplotlib.pyplot as plt

# Grille simple
nlay = 1
nrow = 100
ncol = 100
delr = delc = 10.0
top = 30.0
botm = 0.0
    
def modflow():
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
    wel_spd = {0: [[0, 50, 50, -500]]}  # débit négatif = extraction (m3/jour)
    
    wel = flopy.modflow.ModflowWel(mf, stress_period_data=wel_spd)
    
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
    # plt.figure(figsize=(8, 6))
    # plt.imshow(head[0], cmap="Blues", extent=(0, ncol * delr, 0, nrow * delc))
    # plt.colorbar(label="Niveau piézométrique (m)")
    # plt.scatter((50) * delr, (50) * delc, color="red", label="Zone de perte", s=10)
    # plt.title("Effet d'une perte d'eau locale (puits fictif)")
    # plt.xlabel("Distance Est (m)")
    # plt.ylabel("Distance Nord (m)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    
    return head


def compute_gravity_effect(head, gravimeter_position, cell_size, botm):
    G = 6.67430e-11  # Constante gravitationnelle
    rho = 1000       # Densité de l'eau (kg/m³)

    nlay, nrow, ncol = head.shape
    total_gravity = 0.0

    for row in range(nrow):
        for col in range(ncol):
            # Masse d'eau dans la cellule
            height = max(0, head[0, row, col] - botm)
            volume = height * cell_size ** 2
            mass = rho * volume
 
            if mass == 0:
                continue  # aucune eau, pas de contribution

            # Centre de masse de l'eau
            x = (col + 0.5) * cell_size
            y = (row + 0.5) * cell_size
            z = botm + height / 2  # milieu vertical de l'eau

            # Distance 3D capteur -> centre de masse
            dx = x - gravimeter_position[0]
            dy = y - gravimeter_position[1]
            dz = z - gravimeter_position[2]

            r = np.sqrt(dx**2 + dy**2 + dz**2)

            if r == 0:
                continue  # éviter division par zéro

            # Contribution gravitationnelle
            g = G * mass / r**2
            total_gravity += g

    return total_gravity

def simulate_classical_gravimeter(true_g):
    noise = np.random.normal(0, 5e-8)  # bruit typique ~5 µGal
    return true_g + noise

def simulate_quantum_gravimeter(true_g):
    noise = np.random.normal(0, 1e-8)  # bruit typique ~1 µGal
    bias = np.random.uniform(-2e-8, 2e-8)  # biais laser, calibration
    return true_g + noise + bias


def iteration():
    
    c = 0 
    
    positions = []
    gaps = np.zeros((nrow,ncol,10))
    
    # mat_gap = np.zeros((nrow,ncol))   
    for i in range (2):
        
        positions = []
        
        for row in range(0,nrow,10):
            for col in range(0,ncol,10):
                
                gravimeter_position = [col,row,0]
                g = compute_gravity_effect(modflow(), gravimeter_position, delr, botm)
                
                
                print(str(c*100/((100/5)**2)) + " %")
                
                c += 1
    
                classical_measure = simulate_classical_gravimeter(g)
                quantum_measure = simulate_quantum_gravimeter(g)
                
                classical_error = np.abs(g - classical_measure)
                quantum_error = np.abs(g - quantum_measure)
                
                gap = np.abs(classical_error - quantum_error)
                
                # print(f"classical gravimeter : {classical_measure:.8f} m/s²")
                # print(f"quantum gravimeter : {quantum_measure:.8f} m/s²")
                # print(f"classical gravimeter error : {classical_error:.8f} m/s²")
                # print(f"quantum gravimeter error : {quantum_error:.8f} m/s²")
                # print(f"gap between quantum and classical measurements : {gap:.8f} m/s²")
                
                positions.append([row,col])
                gaps[row,col,i] = gap
                # mat_gap[col,row] = gap
            
    x = [p[0] * delr for p in positions]
    y = [p[1] * delc for p in positions]
    
    gap_mean = np.mean(gaps, axis = 0)
    
    # Visualisation
    ##Groundwater map
    plt.figure(figsize=(8, 6))
    plt.imshow(modflow()[0], cmap="Blues", extent=(0, ncol * delr, 0, nrow * delc))
    plt.colorbar(label="Niveau piézométrique (m)")
    
    ##Gaps points
    # plt.imshow(mat_gap, cmap="viridis", extent=(0, ncol * delr, 0, nrow * delc))
    plt.scatter(x,y , c = gap_mean, cmap = "hot", s=15)
    plt.colorbar(label="gap between quantum and classical (m/s²)")
    # plt.scatter((50) * delr, (50) * delc, color="red", label="Zone de perte", s=10)
    plt.title("Local groundwater loss effect on gravimetric measurements")
    plt.xlabel("East distance (m)")
    plt.ylabel("North distance (m)")
    
    
    plt.grid(False)
    plt.show()
    
    return gaps

if __name__ == "__main__":
    
    
    # gravimeter_position = [ncol * delr / 2, nrow * delc / 2, 0]  # x, y, z
    # gravity = compute_gravity_effect(modflow(), gravimeter_position, delr, botm)

    # print(f"simulated measured gravity : {gravity:.10e} m/s²")
    
    # true_g = gravity + 9.81

    # classical_measure = simulate_classical_gravimeter(true_g)
    # quantum_measure = simulate_quantum_gravimeter(true_g)
    
    # classical_error = np.abs(true_g - classical_measure)
    # quantum_error = np.abs(true_g - quantum_measure)
    
    # print(f"real gravity : {true_g:.8f} m/s²")
    # print(f"classical gravimeter : {classical_measure:.8f} m/s²")
    # print(f"quantum gravimeter : {quantum_measure:.8f} m/s²")
    # print(f"classical gravimeter error : {classical_error:.8f} m/s²")
    # print(f"quantum gravimeter error : {quantum_error:.8f} m/s²")
    
    gaps = iteration()
    