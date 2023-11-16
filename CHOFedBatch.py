import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

t_min, t_max, t_num = 0, 14, 15
step0, step1 = 0, 14
def get_batch(csv_filename):
    batches = []
    batch = []
    with open(csv_filename, newline='') as f:
        reader = csv.reader(f, delimiter=',')
        i = 0
        for k, row in enumerate(reader):
            if k > 0:
                batch.append(np.array(row[1:]).astype('float'))
                i += 1
                if i == 15:
                    batches.append(batch)
                    batch = []
                    i = 0
    batches = np.array(batches)
    return batches

batches = get_batch('in_silico_experiments.csv')

# Constants
mu_max, k_d = 0.029, 0.016  # h^-1
Y_Xv_glc, Y_Xv_gln = 1.69e8, 9.74e8  # cell/mmol
Y_lac_glc, Y_amm_gln = 1.23, 0.67  # mmol/mmol

# Parameters Initial guess
K_glc, K_gln = 0.084, 0.047  # mM
KI_lac, KI_amm = 43, 6.51  # mM
KD_lac, KD_amm = 45.8, 6.51  # mM
m_glc = 69.2e-12  # mmol/cell·h
a1 = 3.2e-12  # mmol/cell·h
a2 = 2.1  # mM
d_gln = 7.2e-3  # h^-1
r_amm = 6.3e-12  # mmol/cell·h
Qp = 4e-12  # L/cell·h

Xv_unit_conv = 10**9

init_params = [K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp]

def dudt(t, U, mu_max, k_d, Y_Xv_glc, Y_Xv_gln, Y_lac_glc, Y_amm_gln, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp, in_glcs, in_glns):
    # Variables
    Xv, GLC, GLN, AMM, LAC, P = U
    # Xv [cell/L], GLC [mM], GLN [mM], LAC [mM], AMM [mM], P [mg/L]
    Xv *= Xv_unit_conv

    in_glc = in_glcs[int(t)]
    in_gln = in_glns[int(t)]

    mu = mu_max * (GLC / (K_glc + GLC)) * (GLN / (K_gln + GLN)) * (KI_lac / (KI_lac + LAC)) * (KI_amm / (KI_amm + AMM))
    mu_d = k_d * (LAC / (KD_lac + LAC)) * (AMM / (KD_amm + AMM))
    m_gln = a1 * GLN / (a2 + GLN)

    dXvdt = (mu - mu_d) * Xv
    dGLCdt = -((mu - mu_d) / Y_Xv_glc + m_glc) * Xv + in_glc
    dGLNdt = -((mu - mu_d) / Y_Xv_gln + m_gln) * Xv - d_gln * GLN + in_gln
    dLACdt = Y_lac_glc * ((mu - mu_d) / Y_Xv_glc + m_glc) * Xv
    dAMMdt = Y_amm_gln * ((mu - mu_d) / Y_Xv_gln) * Xv - r_amm * Xv + d_gln * GLN
    dPdt = Qp * Xv * (1 - mu / mu_max) * P

    return [dXvdt, dGLCdt, dGLNdt, dAMMdt, dLACdt, dPdt]

def mechanistical_model(batch, mu_max, k_d, Y_Xv_glc, Y_Xv_gln, Y_lac_glc, Y_amm_gln, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp):
        U0 = batch[0][1:7]
        tgrid = np.linspace(t_min, t_max, t_num).astype(int)
        in_glcs = np.array(batch[:,9])
        in_glns = np.array(batch[:,10])
        SolU = integrate.solve_ivp(dudt, [step0, step1], U0, t_eval=tgrid, args=(mu_max, k_d, Y_Xv_glc, Y_Xv_gln, Y_lac_glc, Y_amm_gln, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp, in_glcs, in_glns))
        solution = np.array(SolU.y.T)[1:,:].T.flatten()
        return solution

Y = mechanistical_model(batches[15], mu_max, k_d, Y_Xv_glc, Y_Xv_gln, Y_lac_glc, Y_amm_gln, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp)
print('Y:', Y, 'Y.shape', Y.shape)




    # def get_variables(self, Xv, GLC, GLN, LAC, AMM, P):
    #     Xv, GLC, GLN, LAC, AMM, P = self.Xv, self.GLC, self.GLN, self.AMM, self.P
    #     return [Xv, GLC, GLN, LAC, AMM, P]
    #
    #
    # # define objective function
    # def objective_function(self, n):
    #     batch = self.get_batch(n)



    # def mechanistical_model(batch, mu_max, k_d, Y_x_glc, Y_x_gln, Y_lac_glc):
    #     U0 = batch[0][1:7]
    #     tgrid = np.linspace(t_min, t_max, t_num).astype(int)
    #     in_glcs = np.array(batch[:,9])
    #     in_glns = np.array(batch[:,10])
    #     SolU = integrate.solve_ivp(dudt, [step0, step1], U0, t_eval=tgrid, args=(mu_max, k_d, Y_x_glc, Y_x_gln, Y_lac_glc))
    #     solution = np.array(SolU.y.T)[1:,:].T.flatten()
    #     return solution

    # def plotting(self):
    #     # Scetch
    #     mpl.rcParams['lines.markersize'] = 5
    #
    #     # Xv
    #     Xv_cells_mL = SolU.y[0] * 10**-9
    #     plt.scatter(np.array(time_stamps_list), np.array(Xv_list), marker='s', label='Data')
    #     plt.plot(tgrid, Xv_cells_mL, 'ro', label='Model')
    #     plt.xlabel('Time (h)'); plt.ylabel(r'$X_{V}\ (10^{6}\ cells/mL)$')
    #     plt.legend()
    #     plt.show()
    #
    #     # Glucose
    #     plt.scatter(np.array(time_stamps_list), np.array(GLC_list), marker='s', label='Data')
    #     plt.plot(tgrid, SolU.y[1], 'ro', label='Model')
    #     plt.xlabel('Time (h)'); plt.ylabel('Glucose (mM)')
    #     plt.legend()
    #     plt.show()
    #
    #     # Glutamine
    #     plt.scatter(np.array(time_stamps_list), np.array(GLN_list), marker='s', label='Data')
    #     plt.plot(tgrid, SolU.y[2], 'ro', label='Model')
    #     plt.xlabel('Time (h)'); plt.ylabel('Glutamine (mM)')
    #     plt.legend()
    #     plt.show()
    #
    #     # Lactate
    #     plt.scatter(np.array(time_stamps_list), np.array(LAC_list), marker='s', label='Data')
    #     plt.plot(tgrid, SolU.y[3], 'ro', label='Model')
    #     plt.xlabel('Time (h)'); plt.ylabel('Lactate (mM)')
    #     plt.legend()
    #     plt.show()
    #
    #     # Ammonia
    #     plt.scatter(np.array(time_stamps_list), np.array(AMM_list), marker='s', label='Data')
    #     plt.plot(tgrid, SolU.y[4], 'ro', label='Model')
    #     plt.xlabel('Time (h)');
    #     plt.ylabel('Ammonia (mM)')
    #     plt.legend()
    #     plt.show()
    #
    #     # Product
    #     plt.scatter(np.array(time_stamps_list), np.array(P_list), marker='s', label='Data')
    #     plt.plot(tgrid, SolU.y[5], 'ro', label='Model')
    #     plt.xlabel('Time (h)'); plt.ylabel('Product (mg/L)')
    #     plt.legend()
    #     plt.show()
