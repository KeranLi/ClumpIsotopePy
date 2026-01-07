import numpy as np
import isotopylog as ipl

MYR_TO_SEC = 1e6 * 365 * 24 * 60 * 60


def compute_history(time_myr, T_k, ed, d0_std):
    time_myr = np.asarray(time_myr, dtype=float)
    T_k = np.asarray(T_k, dtype=float)

    time_rev = time_myr[::-1]
    T_rev = T_k[::-1]

    t_sec = time_rev * MYR_TO_SEC
    t_sec = t_sec[0] - t_sec

    D0 = ipl.Deq_from_T(T_rev[0])
    d0 = [D0, 0, 0]

    D_rev, Dstd_rev = ipl.geologic_history(t_sec, T_rev, ed, d0, d0_std=[d0_std])

    D = D_rev[::-1]
    Dstd = Dstd_rev[::-1]
    Deq = ipl.Deq_from_T(T_k)

    return D, Dstd, Deq
