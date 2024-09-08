import numpy as np
import matplotlib.pyplot as plt

schemas = {1: 'Centré', 2: 'Décentré amont', 3: 'Lax-Friedrichs', 4: 'Lax-Wendroff'}

# Fonction pour la solution exacte
def solution_burgers(x, t):
    return np.where(((x-2)*(1/t) < 0.25), 0.4, 0.1)
    
    
# Fonction pour la solution initiale
def solution_initiale(x):
    return np.where((x < 2), 0.4, 0.1)



def simulation_burgers(num_schema, L=10,N=100,CFL=0.8,tfinal=6):

    # Paramètres
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    u = np.zeros(N)
    unew = np.zeros(N)
    
    # Fonction pour la solution initiale
    u = solution_initiale(x)

    # Initialisation des listes pour stocker les résultats
    Ts = []
    Us = []
    

    temps = 0
    
    while temps <= tfinal:
        max_u = max(u)
        dt = CFL * dx / max_u
        # lamda = max_u * dt / dx
        if num_schema == 1:
            for i in range(1, N - 1):
                unew[i] = u[i] - 0.5 * (u[i] * dt / dx) * (u[i + 1] - u[i - 1])
        elif num_schema == 2:
            for i in range(1, N - 1):
                unew[i] = u[i] - (u[i] * dt / dx) * (u[i] - u[i - 1])
        elif num_schema == 3:
            for i in range(1, N - 1):
                unew[i] = 0.5 * (u[i + 1] + u[i - 1]) - 0.5*(dt / (2 * dx)) *(u[i + 1] + u[i - 1])* (u[i + 1] - u[i - 1])
        elif num_schema == 4:
            # Lax-Wendroff
            for i in range(1, N - 1):
                unew[i] = u[i] - (u[i]**2)*(dt / (2 * dx)) * (u[i + 1] - u[i - 1]) + (u[i]**2)*(dt**2 / (2 * dx**2)) * (u[i + 1] - 2 * u[i] + u[i - 1])
                
        unew[0] = unew[1]
        unew[N - 1] = unew[N - 2]
        u = unew.copy()
        Ts.append(temps)
        Us.append(u)
        temps += dt
    

    # intervalle autour de 2.4 et 4.4 pour l'affichage
    t_2_5 = max([t for t in Ts if 2.4 <= t <= 2.5])
    t_4_5 = max([t for t in Ts if 4 <= t <= 4.5])
    # Tracé des résultats
    plt.figure(figsize=(12, 6))
    plt.suptitle(f'Comparaison des solutions numériques et exactes pour {schemas[num_schema]}')
    plt.subplot(2, 2, 1)
    plt.plot(x, solution_burgers(x, t_2_5), label='Solution exacte', color='blue')
    plt.plot(x, Us[Ts.index(t_2_5)], label=f'Solution numérique t = {t_2_5:.2f}s', linestyle='dashed', color='red')
    plt.title(f'Solution pour N={N} (2.4 <= t <= 2.5)')
    plt.legend()
    plt.grid()
    plt.ylim(0, 1)

    plt.subplot(2, 2, 2)
    plt.plot(x, solution_burgers(x, t_4_5), label='Solution exacte', color='blue')
    plt.plot(x, Us[Ts.index(t_4_5)], label=f'Solution numérique t = {np.round(t_4_5, 3)}s', linestyle='dashed', color='red')
    plt.title(f'Solution pour N={N} (4.4 <= t <= 4.5)')
    plt.legend()
    plt.grid()
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()
    err_2_5 = np.sum(np.abs(solution_burgers(x, t_2_5)-Us[Ts.index(t_2_5)]))
    err_4_5 = np.sum(np.abs(solution_burgers(x, t_4_5)-Us[Ts.index(t_4_5)]))
    print(f"Pour t=2.5, erreur = {err_2_5}")
    print(f"Pour t=4.5, erreur = {err_4_5}")
    
    return


#simulation_burgers(1, N=100)

simulation_burgers(2, N=100)

#simulation_burgers(3, N=100)

#simulation_burgers(4, N=100)