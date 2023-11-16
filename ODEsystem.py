# def dudt(t, U):
#     # Parameters
#     mu_max, k_d = 0.029, 0.016  # h^-1
#     Y_Xv_glc, Y_Xv_gln = 1.69e8, 9.74e8  # cell/mmol
#     Y_lac_glc, Y_amm_gln = 1.23, 0.67  # mmol/mmol
#     K_glc, K_gln = 0.084, 0.047  # mM
#     KI_lac, KI_amm = 43, 6.51  # mM
#     KD_lac, KD_amm = 45.8, 6.51  # mM
#     m_glc = 69.2e-12  # mmol/cell·h
#     a1 = 3.2e-12  # mmol/cell·h
#     a2 = 2.1  # mM
#     d_gln = 7.2e-3  # h^-1
#     r_amm = 6.3e-12  # mmol/cell·h
#     Qp = 4e-12  # L/cell·h
#     c_glc, c_gln = 0, 0  # mM
#
#     # Variables
#     Xv, GLC, GLN, LAC, AMM, P = U
#     # Xv [cell/L], GLC [mM], GLN [mM], LAC [mM], AMM [mM], P [mg/L]
#
#     if np.round(t) == 0:
#         c_glc = Feed_GLC_list[0]
#         c_gln = Feed_GLN_list[0]
#
#     if t > 1:
#         for time in time_stamps_list:
#             if t_cache[-1] <= time and time <= t:
#                 print('previous time:', t_cache[-1], 'current time:', t)
#                 index = time_stamps_list.index(time)
#                 c_glc = Feed_GLC_list[index]
#                 c_gln = Feed_GLN_list[index]
#
#     mu = mu_max * (GLC / (K_glc + GLC)) * (GLN / (K_gln + GLN)) * (KI_lac / (KI_lac + LAC)) * (KI_amm / (KI_amm + AMM))
#     mu_d = k_d * (LAC / (KD_lac + LAC)) * (AMM / (KD_amm + AMM))
#
#     if AMM <= 4.6:
#         mu = mu_max * (GLC / (K_glc + GLC)) * (GLN / (K_gln + GLN)) * (KI_lac / (KI_lac + LAC))
#
#     if LAC <= 52:
#         mu = mu_max * (GLC / (K_glc + GLC)) * (GLN / (K_gln + GLN)) * (KI_amm / (KI_amm + AMM))
#
#     m_gln = a1 * GLN / (a2 + GLN)
#
#     if t <= 144:
#         m_glc, m_gln = 0, 0
#
#     dXvdt = (mu - mu_d) * Xv
#     dGLCdt = -((mu - mu_d) / Y_Xv_glc + m_glc) * Xv + c_glc
#     dGLNdt = -((mu - mu_d) / Y_Xv_gln + m_gln) * Xv - d_gln * GLN + c_gln
#     dLACdt = Y_lac_glc * ((mu - mu_d) / Y_Xv_glc + m_glc) * Xv
#     dAMMdt = Y_amm_gln * ((mu - mu_d) / Y_Xv_gln) * Xv - r_amm * Xv + d_gln * GLN
#     dPdt = Qp * Xv * (1 - mu / mu_max) * P
#
#     t_cache.append(t)
#
#     if Xv < 0:
#         Xv = 0
#     if GLC < 0:
#         GLC = 0
#     if GLN < 0:
#         GLN = 0
#     if LAC < 0:
#         LAC = 0
#     if AMM < 0:
#         AMM = 0
#     if P < 0:
#         P = 0
#
#     return [dXvdt, dGLCdt, dGLNdt, dLACdt, dAMMdt, dPdt]
#
# U0 = [float(batch[0][2])*10**9, float(batch[0][3]), float(batch[0][4]), float(batch[0][5]), float(batch[0][6]), float(batch[0][7])]
#         SolU = integrate.solve_ivp(dudt, [t_min, t_max], U0, method='RK23', t_eval=tgrid)