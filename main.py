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
in_glc, in_gln = 0, 0  # mM

Xv_unit_conv = 10**9

init_params = [K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp]

def dudt(t, U, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp, in_glcs, in_glns):
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

def mechanistical_model(batch, K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp):
        U0 = batch[0][1:7]
        tgrid = np.linspace(t_min, t_max, t_num).astype(int)
        in_glcs = np.array(batch[:,9])
        in_glns = np.array(batch[:,10])
        SolU = integrate.solve_ivp(dudt, [step0, step1], U0, t_eval=tgrid, args=[K_glc, K_gln, KI_lac, KI_amm, KD_lac, KD_amm, m_glc, a1, a2, d_gln, r_amm, Qp, in_glcs, in_glns])
        solution = np.array(SolU.y.T)[1:,:].T.flatten()
        return solution

popt = []
for k in range(78):
    try:
        if k == 0:
            popt, pcov = curve_fit(mechanistical_model, batches[k], batches[k][1:, 1:7].T.flatten(), p0=init_params)
            popt = popt[:,np.newaxis]
        else:
            popt_other, pcov = curve_fit(mechanistical_model, batches[k], batches[k][1:, 1:7].T.flatten(), p0=init_params,\
                                   bounds=([0.01, 0.01, 30, 5, 40, 5, 50e-12 * 24, 1e-12 * 24, 1, 5e-3 * 24, 0.01e-12 * 24, 0.001e-12 * 24],
                                           [0.2, 0.2, 70, 20, 70, 20, 90e-12 * 24, 10e-12 * 24, 10, 20e-3 * 24, 20e-12 * 24, 10e-12 * 24]))
            popt = np.concatenate((popt, popt_other[:,np.newaxis]), axis=1)
        print(f'{k+1}th :', popt[:,-1].T)
    except:
        continue
popt = np.mean(popt)

print('final params values:', popt)


def get_mse(diff):
    return np.mean(np.square(diff), axis=1)

for k in range(78, 86):
    labels = (batches[k][1:][1:7]).T
    preds = mechanistical_model(batches[k], *popt).reshape(5, 14)
    difference = preds - labels
    if k == 78:
        diff_list = difference
    else:
        diff_list = np.concatenate((diff_list, difference), axis=1)

mse = get_mse(diff_list)
rmse = np.sqrt(mse)
order = ['Xv', 'GLC', 'GLN', 'AMM', 'LAC', 'P']
print('total loss :', np.sum(mse))

print('<MSE>')
for l in range(len(mse)):
    print(order[l], ':', mse[l])

print('<RMSE>')
for l in range(len(rmse)):
    print(order[l], ':', rmse[l])












