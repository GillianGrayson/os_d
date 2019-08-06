#include "characteristics.h"

void characteristics_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id)
{
	// Add here regular characteristics
}
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id)
{
	MKL_Complex16 * rho_in_d = new MKL_Complex16[md.size * md.size];
	MKL_Complex16 * evals = new MKL_Complex16[md.size];

	for (int state_id_1 = 0; state_id_1 < md.size; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < md.size; state_id_2++)
		{
			rho_in_d[state_id_1 * md.size + state_id_2].real = 0.0;
			rho_in_d[state_id_1 * md.size + state_id_2].imag = 0.0;
		}
	}

	for (int i = 0; i < m->Rho->N; i++)
	{
		for (int k = m->Rho->RowIndex[i]; k < m->Rho->RowIndex[i + 1]; k++)
		{
			rho_in_d[i * md.size + m->Rho->Col[k]].real = m->Rho->Value[k].re;
			rho_in_d[i * md.size + m->Rho->Col[k]].imag = m->Rho->Value[k].im;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md.size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md.size; st_id_2++)
		{
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].real += rho_in_d[st_id_1 * md.size + st_id_2].real;
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].imag += rho_in_d[st_id_1 * md.size + st_id_2].imag;
		}
	}

	int info = LAPACKE_zgeev(
		LAPACK_ROW_MAJOR,
		'N',
		'N',
		md.size,
		rho_in_d,
		md.size,
		evals,
		NULL,
		md.size,
		NULL,
		md.size);

	if (info > 0)
	{
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}

	for (int st_id = 0; st_id < md.size; st_id++)
	{
		pd.deep_evals[dump_id * md.size + st_id] = evals[st_id].real;
	}

	delete[] evals;
	delete[] rho_in_d;
}