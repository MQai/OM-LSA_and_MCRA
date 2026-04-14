#include "OMLSA-MCRA.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define min(a,b)            (((a) < (b)) ? (a) : (b))
using namespace std;

[[noreturn]] static void mcra_fatal(const char* msg) {
	cerr << "MCRA fatal: " << msg << endl;
	abort();
}

Cordic::Cordic()
{
}
int Cordic::cordic_sin(int theta){
	int negtive_flag = 1; 
	int sin_value, d;
	x[0] = 652032943;
	y[0] = 0;
	while (theta < -Pi_2) {
		theta += Pi_2; 
		negtive_flag *= -1; 
	}
	while (theta > Pi_2) {
		theta -= Pi_2; 
		negtive_flag *= -1; 
	}
	while (theta > Pi_div2) {
		theta = Pi_2 - theta; 
	}
	while (theta < -Pi_div2) {
		theta = -Pi_2 - theta; 
	}

	for (int i = 0; i < 15; i++) {
		if (theta > 0)d = 1;
		else d = -1;
		x[i + 1] = x[i] - d * (y[i] >> i);
		y[i + 1] = y[i] + d * (x[i] >> i);
		theta = theta - d * atanget[i];
	}
	sin_value = negtive_flag * y[15];
	return sin_value;
}

int Cordic::cordic_cos(int theta)
{
	int negtive_flag = 1;
	int cos_value, d;
	x[0] = 652032943;
	y[0] = 0;
	while (theta < -Pi_2) {
		theta += Pi_2;
		negtive_flag *= -1;
	}
	while (theta > Pi_2) {
		theta -= Pi_2;
		negtive_flag *= -1;
	}
	while (theta > Pi_div2) {
		theta = Pi_2 - theta;
		negtive_flag *= -1;
	}
	while (theta < -Pi_div2) {
		theta = -Pi_2 - theta;
		negtive_flag *= -1;
	}

	for (int i = 0; i < 15; i++) {
		if (theta > 0)d = 1;
		else d = -1;
		x[i + 1] = x[i] - d * (y[i] >> i);
		y[i + 1] = y[i] + d * (x[i] >> i);
		theta = theta - d * atanget[i];
	}
	cos_value = negtive_flag * x[15];
	return cos_value;
}

Cordic::~Cordic()
{
}

MY_B4_FFT::MY_B4_FFT() {
}

void MY_B4_FFT::initial(int N) {
	m_fwlen = N;
	m_finc = m_fwlen>>1;
	m_M4 = log(m_fwlen) / log(4);
	int N_pow = pow(4, m_M4);
	if (N_pow - m_fwlen!=0  || m_fwlen<=0|| m_fwlen> 16384) {
		mcra_fatal("MY_B4_FFT::initial: N must be a power of 4 in (0, 16384]");
	}
	m_Buffer_cos = new int[m_fwlen];
	m_Buffer_sin = new int[m_fwlen];
	m_sort4_count = new int[m_fwlen];
	m_sort_temp_r = new int[m_fwlen];
	m_sort_temp_i = new int[m_fwlen];
	m_e1_position = new int[m_fwlen];
	m_o1_position = new int[m_fwlen];
	m_e2_position = new int[m_fwlen];
	m_o2_position = new int[m_fwlen];
	m_yr = new int[m_fwlen];
	m_yi = new int[m_fwlen];
	m_ifft_move_bit = m_M4<<1;

	for (int i = 0; i < m_fwlen; i++) {  
		m_Buffer_cos[i] = (coc.cordic_cos((int)(-2 * Pi*i * 32768)>>m_ifft_move_bit));
		m_Buffer_sin[i] = (coc.cordic_sin((int)(-2 * Pi*i * 32768)>>m_ifft_move_bit));
	}
	Base4_Sort();
}

void MY_B4_FFT::Base4_Sort(){
	int  bit_rest, space, add;
	memset(m_sort4_count, 0, sizeof(int)*m_fwlen);
	for (int l = 0; l < m_M4; l++)
	{
		space = pow(4, l);
		add = pow(4, m_M4 - l - 1);
		for (int i = 0; i < m_fwlen; i++){
			bit_rest = (i / space) % 4;
			if (bit_rest != 0)			   
				m_sort4_count[i] += add * bit_rest;
		}
	}
}

void MY_B4_FFT::base4_fft(Complex_num *x, int sign) {
	int i, j, k, tr, ti, wi2, wi3;
	if (x == NULL) {
		mcra_fatal("MY_B4_FFT::base4_fft: null input buffer");
	}
	if (!abs(sign)) {
		mcra_fatal("MY_B4_FFT::base4_fft: sign must be non-zero");
	}
	int ctemp;
	for (int i = 0; i < m_fwlen; i++){
		m_sort_temp_r[i] = x[i].real;
		m_sort_temp_i[i] = x[i].imag;
	}
	for (int i = 0; i < m_fwlen; i++){
		ctemp = m_sort4_count[i];
		x[i].real = m_sort_temp_r[ctemp];
		x[i].imag = m_sort_temp_i[ctemp];
	}

	int BlockLen, BlockNum, BlockLen2;  int wi, F1, F2; int F3, F0; Complex_num X0, X1, X2, X3;
	for (i = 1; i <= m_M4; i++) {
		BlockNum = pow(4, m_M4 - i);
		BlockLen = pow(4, i);
		BlockLen2 = BlockLen >> 2;
		for (j = 0; j < BlockNum; j++) {
			for (k = 0; k < BlockLen2; k++) {
				F0 = k + j * BlockLen;
				F1 = F0 + BlockLen2;
				F2 = F1 + BlockLen2;
				F3 = F2 + BlockLen2;
				wi = BlockNum * k;
				wi2 = wi << 1;
				wi3 = (wi << 1) + wi;

				X0.real = x[F0].real;
				X0.imag = x[F0].imag;
				X1.real = (((int64_t)x[F1].real*m_Buffer_cos[wi]) >> 30) - (((int64_t)x[F1].imag*m_Buffer_sin[wi]) >> 30);
				X1.imag = (((int64_t)x[F1].imag*m_Buffer_cos[wi]) >> 30) + (((int64_t)x[F1].real*m_Buffer_sin[wi]) >> 30);
				X2.real = (((int64_t)x[F2].real*m_Buffer_cos[wi2]) >> 30) - (((int64_t)x[F2].imag*m_Buffer_sin[wi2]) >> 30);
				X2.imag = (((int64_t)x[F2].imag*m_Buffer_cos[wi2]) >> 30) + (((int64_t)x[F2].real*m_Buffer_sin[wi2]) >> 30);
				X3.real = (((int64_t)x[F3].real*m_Buffer_cos[wi3]) >> 30) - (((int64_t)x[F3].imag*m_Buffer_sin[wi3]) >> 30);
				X3.imag = (((int64_t)x[F3].imag*m_Buffer_cos[wi3]) >> 30) + (((int64_t)x[F3].real*m_Buffer_sin[wi3]) >> 30);

				x[F0].real = X0.real + X1.real + X2.real + X3.real;
				x[F0].imag = X0.imag + X1.imag + X2.imag + X3.imag;
				x[F1].real = X0.real + X1.imag - X2.real - X3.imag;
				x[F1].imag = X0.imag - X1.real - X2.imag + X3.real;
				x[F2].real = X0.real - X1.real + X2.real - X3.real;
				x[F2].imag = X0.imag - X1.imag + X2.imag - X3.imag;
				x[F3].real = X0.real - X1.imag - X2.real + X3.imag;
				x[F3].imag = X0.imag + X1.real - X2.imag - X3.real;
			}
		}
	}

	if (sign == -1)
	{
		x[0].real = x[0].real >> m_ifft_move_bit;
		x[0].imag = x[0].imag >> m_ifft_move_bit;
		for (i = 1; i <= m_finc; i++) {
			m_fly_tempr = (x[m_fwlen - i].real) >> m_ifft_move_bit;
			m_fly_tempi = (x[m_fwlen - i].imag) >> m_ifft_move_bit;
			x[m_fwlen - i].real = (x[i].real) >> m_ifft_move_bit;
			x[m_fwlen - i].imag = (x[i].imag) >> m_ifft_move_bit;
			x[i].real = m_fly_tempr;
			x[i].imag = m_fly_tempi;
		}
	}
	else if (sign == 1) {  
		for (int i = 0; i < m_fwlen; i++) {
			if (i == 0) {
				m_yr[i] = x[i].real; 
				m_yi[i] = x[i].imag;
			}
			else {
				m_yr[i] = x[m_fwlen - i].real; 
				m_yi[i] = x[m_fwlen - i].imag;
			}
			m_e1_position[i] = (x[i].real + m_yr[i]) >> 1;
			m_o1_position[i] = (x[i].real - m_yr[i]) >> 1;
			m_e2_position[i] = (x[i].imag + m_yi[i]) >> 1;
			m_o2_position[i] = (x[i].imag - m_yi[i]) >> 1;
		}
		for (int i = 0; i < m_fwlen; i++) {
			x[i].real = m_e1_position[i];
			x[i].imag = m_o2_position[i];
			x[i + m_fwlen].real = m_e2_position[i];
		 	x[i + m_fwlen].imag = -m_o1_position[i];
		}
	}
}

MY_B4_FFT::~MY_B4_FFT() {
	delete[] m_Buffer_cos;
	m_Buffer_cos = nullptr;
	delete[] m_Buffer_sin;
	m_Buffer_sin = nullptr;
	delete[] m_sort4_count;
	m_sort4_count = nullptr;
	delete[] m_sort_temp_r;
	m_sort_temp_r = nullptr;
	delete[] m_sort_temp_i;
	m_sort_temp_i = nullptr;
	delete[] m_e1_position;
	m_e1_position = nullptr;
	delete[] m_o1_position;
	m_o1_position = nullptr;
	delete[] m_yr;
	m_yr = nullptr;
	delete[] m_e2_position;
	m_e2_position = nullptr;
	delete[] m_o2_position;
	m_o2_position = nullptr;
	delete[] m_yi;
	m_yi = nullptr;
}

void MY_B4_FFT::fftfix(Complex_num *x, int sign) {
	int i, j, k, u = 0, l = 0, wi = 0, n1, tr, ti, N = m_fwlen;
	int SubBlockNum, SubBlockStep = 1;

	int tempr, tempi;
	N = N << 1;

	n1 = m_fwlen - 1;
	for (j = 0, i = 0; i < n1; i++)
	{
		if (i < j)
		{
			tr = x[i].real;
			ti = x[i].imag;

			x[i].real = x[j].real;
			x[i].imag = x[j].imag;

			x[j].real = tr;
			x[j].imag = ti;
		}
		k = N / 2;
		while (k < (j + 1))
		{
			j = j - k;
			k = k / 2;
		}
		j = j + k;
	}

	for (k = N; k > 1; k = (k >> 1)) {
		SubBlockNum = k >> 1;
		SubBlockStep = SubBlockStep << 1;
		wi = 0;
		for (j = 0; j < SubBlockStep >> 1; j++) {
			for (u = j; u < N; u += SubBlockStep) {
				l = u + (SubBlockStep >> 1);
				tempr = (((int64_t)x[l].real*m_Buffer_cos[wi]) >> 30) - (((int64_t)x[l].imag*m_Buffer_sin[wi]) >> 30);
				tempi = (((int64_t)x[l].imag*m_Buffer_cos[wi]) >> 30) + (((int64_t)x[l].real*m_Buffer_sin[wi]) >> 30);

				x[l].real = x[u].real - tempr;
				x[l].imag = x[u].imag - tempi;
				x[u].real = x[u].real + tempr;
				x[u].imag = x[u].imag + tempi;
			}
			wi += SubBlockNum;
		}
	}

	if (sign == -1)
	{
		for (i = 0; i < N; i++)
		{
			x[i].real /= N;
			x[i].imag /= N;
		}
	}
}

G_calculate::G_calculate()
{
}
void G_calculate::Initialize(int wlen) {
	m_lwlen = wlen ;
	m_linc = m_lwlen >> 1;
	m_linc_move = log(m_linc) / log(2);
	N_wlen1 = m_lwlen + 1;
	m_num_mag = 128;
	m_num_mag_pow = 16384;
	m_num_mag_pow2 = 268435456;
	m_amp_para = 7, m_amp_para_double = 14, m_amp_para_three = 21;  
	m_snpramin = 0.03*m_num_mag_pow;
	beta = 0.7*m_num_mag_pow;
	m_cosen_max = 0.8*m_num_mag_pow;
	m_cosen_min = 0.1*m_num_mag_pow;
	m_cosen_pmax = 5 * m_num_mag_pow;
	m_cosen_pmin = 1 * m_num_mag_pow;
	m_cosen_max_min = log(8) *m_num_mag_pow;
	m_w_global = 8; m_wt_global = 17;

	m_arr_temp = new int[m_lwlen];
	m_S_f = new int[m_lwlen];
	m_abs_Y = new unsigned int[m_lwlen];
	m_dataBuf = new short[m_lwlen];
	m_init_S = new int[m_lwlen];
	m_Gh1 = new int[m_lwlen];
	m_q = new int[m_lwlen];
	m_pp = new int[m_lwlen];
	m_G = new int[m_lwlen];
	m_M = new unsigned int[m_lwlen];
	m_ns_storage = new short[m_lwlen];
	m_cosen_local = new int[m_lwlen];
	m_cosen_global = new int[m_lwlen];
	m_plocal = new int[m_lwlen];
	m_pglobal = new int[m_lwlen];
	m_h_global = new short[m_wt_global];
	m_old_cosen = new int[m_lwlen];
	m_cosen = new int[m_lwlen];
	m_lamda_d = new int[N_wlen1];
	m_integra = new int[m_lwlen];
	m_init_S_min = new int[m_lwlen];
	m_init_S_tmp = new int[m_lwlen];
	m_init_p1 = new int[m_lwlen];
	m_EN_cos = new short[m_lwlen];
	m_EN_sin = new short[m_lwlen];
	m_post_SNR = new int[m_lwlen];
	m_pr_SNR = new int[m_lwlen];
	m_E_pr_SNR = new int[m_lwlen];
	m_v = new int[m_lwlen];
	
	for (int i = 1; i < (m_wt_global + 1); ++i)
		m_h_global[i - 1] = (0.5 - 0.5 * cos(2.0 * Pi*(i) / (m_wt_global + 1))) * 16384;

	for (int i = 0; i < m_lwlen; i++) 
		m_pr_SNR[i] = 0.98*m_num_mag_pow;

	memset(m_old_cosen, 0, sizeof(int)*m_lwlen);
	memset(m_cosen, 0, sizeof(int)*m_lwlen);
	memset(m_dataBuf, 0, sizeof(short) * m_lwlen);
	memset(m_ns_storage, 0, sizeof(short)*m_lwlen);

}

void G_calculate::NoiseEstimation(int blockInd)
{
	int p;
	int L;
	if (blockInd >= 10) {
		L = 50;
	}
	else {
		L = 1;
		if (blockInd == 0) {
			for (int i = 0; i <= m_linc; i++) {
				m_init_S[i] = m_abs_Y[i];
				m_init_S_min[i] = m_init_S[i];
				m_init_S_tmp[i] = m_init_S[i];
				m_lamda_d[i] = m_init_S[i];
			}
			memset(m_init_p1, 0, sizeof(int)*(m_linc + 1));
		}
	}

	if ((blockInd + 1) % L == 0) {
		for (int k = 0; k <= m_linc; k++) {
			m_init_S_min[k] = min(m_init_S_tmp[k], m_init_S[k]);
			m_init_S_tmp[k] = m_init_S[k];
		}
	}

	for (int k = 0; k <= m_linc; k++)
	{
		m_init_S[k] = (m_init_S[k] >> 1) + (m_init_S[k] >> 2) + (m_init_S[k] >> 4) + (m_abs_Y[k] >> 3) + (m_abs_Y[k] >> 4);
		m_init_S_min[k] = min(m_init_S_min[k], m_init_S[k]);
		m_init_S_tmp[k] = min(m_init_S_tmp[k], m_init_S[k]);

		if (m_init_S[k] > (m_init_S_min[k] >> 2) + (m_init_S_min[k] << 1))
			p = m_num_mag_pow;
		else
			p = 0;

		m_init_p1[k] = (m_init_p1[k] >> 2) + (p >> 1) + (p >> 2);
		m_lamda_d[k] = ((int64_t)((m_lamda_d[k] >> 1) + (m_lamda_d[k] >> 2) + (m_lamda_d[k] >> 3) + (m_lamda_d[k] >> 4) + (m_abs_Y[k] >> 4))
			*(m_num_mag_pow - m_init_p1[k]) + (int64_t)m_lamda_d[k] * m_init_p1[k]) >> m_amp_para_double;
		m_lamda_d[m_lwlen - k] = m_lamda_d[k];
	}
}
void G_calculate::SpeechAbsenceEstm()
{
	unsigned long long sum = 0, a = 0;
	int p_frame, mu, cosen_peak;
	short old_cosen_frame = 0, cosen_frame = 0;

	for (int k = 0; k <= (m_linc + m_w_global); k++) {
		m_cosen[k] = (m_old_cosen[k] >> 1) + (m_old_cosen[k] >> 2) + (m_E_pr_SNR[k] >> 2);
	}

	for (int k = 0; k <= m_linc; k++) {
		if (k <= m_w_global - 1) {
			m_cosen_global[k] = m_cosen[k];
			if (k == 0)
				m_cosen_local[k] = m_cosen[k];
			else
				m_cosen_local[k] = (m_cosen[k - 1] + (m_cosen[k] << 1) + m_cosen[k + 1]) >> 2;
		}
		else {
			a = 0;
			for (int j = 0; j < 2 * m_w_global + 1; j++)
				a = a + ((int64_t)m_h_global[j] * m_cosen[k + m_w_global - j]);
			m_cosen_global[k] = ((a >> m_amp_para_double) / (m_w_global + 1));
			m_cosen_local[k] = (m_cosen[k - 1] + (m_cosen[k] << 1) + m_cosen[k + 1]) >> 2;
		}

		if (m_cosen_local[k] <= m_cosen_min)
			m_plocal[k] = 0;
		else if (m_cosen_local[k] >= m_cosen_max)
			m_plocal[k] = m_num_mag_pow;
		else
			m_plocal[k] = m_num_mag_pow2 * (log(m_cosen_local[k]) - log(m_cosen_min)) / m_cosen_max_min;

		if (m_cosen_global[k] <= m_cosen_min)
			m_pglobal[k] = 0;
		else if (m_cosen_global[k] >= m_cosen_max)
			m_pglobal[k] = m_num_mag_pow;
		else
			m_pglobal[k] = m_num_mag_pow2 * (log(m_cosen_global[k]) - log(m_cosen_min)) / m_cosen_max_min;
		sum += m_cosen[k];
	}

	cosen_frame = sum >> m_linc_move;
	sum = 0;
	cosen_peak = min(max(cosen_frame, m_cosen_pmin), m_cosen_pmax);

	if (cosen_frame <= ((cosen_peak * m_cosen_min) >> m_amp_para_double))
		mu = 0;
	else if (cosen_frame >= ((cosen_peak * m_cosen_max) >> m_amp_para_double))
		mu = m_num_mag_pow;
	else
		mu = m_num_mag_pow * log(cosen_frame  * m_num_mag_pow / cosen_peak / m_cosen_min) / log(m_cosen_max / m_cosen_min);

	if (cosen_frame > m_cosen_min)
	{
		if (cosen_frame > old_cosen_frame)
			p_frame = m_num_mag_pow;
		else
			p_frame = mu;
		//cout<<"p_frame="<<p_frame<<endl;
	}
	else
	{
		p_frame = 0;
	}
	for (int k = 0; k <= m_linc; k++)
	{
		m_q[k] = m_num_mag_pow - ((int64_t)m_plocal[k] * m_pglobal[k] * p_frame >> 28);
		m_q[k] = min(m_q[k], 0.95*m_num_mag_pow);
		m_old_cosen[k] = m_cosen[k];
		m_old_cosen[m_lwlen - k] = m_old_cosen[k];
	}
}
void G_calculate::G_calculate_process(Complex_num* winData, int blockInd) {
	if (winData == NULL)
		mcra_fatal("G_calculate::G_calculate_process: null winData");

	int post_temp, w = 8;
	for (int i = 0; i <= m_linc + w; i++) {
		m_abs_Y[i] = sqrt(pow(winData[i].real, 2) + pow(winData[i].imag, 2));
		m_EN_cos[i] = ((int64_t)winData[i].real << 12) / (1 > m_abs_Y[i] ? 1 : m_abs_Y[i]);
		m_EN_sin[i] = ((int64_t)winData[i].imag << 12) / (1 > m_abs_Y[i] ? 1 : m_abs_Y[i]);
	}
	NoiseEstimation(blockInd);
	for (int i = 0; i <= m_linc; i++) {
		m_post_SNR[i] = min((int64_t)pow((((int64_t)m_abs_Y[i] * m_num_mag_pow) / (1 > m_lamda_d[i] ? 1 : m_lamda_d[i])), 2) >> m_amp_para_double, 4096 << 14);
		post_temp = max(m_post_SNR[i] - m_num_mag_pow, 0);
		m_E_pr_SNR[i] = min(max((m_pr_SNR[i] >> 1) + (m_pr_SNR[i] >> 2) + (m_pr_SNR[i] >> 3) + (post_temp >> 3), m_snpramin), 4096 << m_amp_para_double);
		m_E_pr_SNR[m_lwlen - i] = m_E_pr_SNR[i];
		m_v[i] = max(min((int64_t)((int64_t)m_E_pr_SNR[i] * m_post_SNR[i] << 10) / (m_num_mag_pow + m_E_pr_SNR[i]), 15 << 24), 1678);
		m_integra[i] = expintpow_solution(m_v[i]);
		m_Gh1[i] = min((int64_t)m_E_pr_SNR[i] * m_integra[i] / (m_num_mag_pow + m_E_pr_SNR[i]), 70 << 14);
	}
	SpeechAbsenceEstm();
	for (int i = 0; i <= m_linc; i++) {
		m_integra[i] = subexp_solution(m_v[i]);
		m_arr_temp[i] = m_num_mag_pow + ((int64_t)(m_num_mag_pow + m_E_pr_SNR[i]) * m_integra[i] * m_q[i] >> m_amp_para_double) / (m_num_mag_pow - m_q[i]);

		m_pp[i] = min(m_num_mag_pow2 / m_arr_temp[i], 1 << 14);
		m_G[i] = Gvalue_solution(m_Gh1[i], m_pp[i]);
		m_M[i] = ((int64_t)m_G[i] * m_abs_Y[i]) >> m_amp_para_double; 
		winData[i].real = ((int64_t)m_M[i] * m_EN_cos[i]) >> 12; 
		winData[i].imag = ((int64_t)m_M[i] * m_EN_sin[i]) >> 12;
		m_pr_SNR[i] = min(pow((int64_t)m_M[i] * m_num_mag / (1 > m_lamda_d[i] ? 1 : m_lamda_d[i]), 2), 4096 << m_amp_para_double);  // 10000
		winData[m_lwlen - i].real = winData[i].real;
		winData[m_lwlen - i].imag = -winData[i].imag;
		aa[i] = m_G[i];
	}
}


G_calculate::~G_calculate()
{
	delete[] m_arr_temp;
	m_arr_temp = nullptr;
	delete[] m_S_f;
	m_S_f = nullptr;
	delete[] m_abs_Y;
	m_abs_Y = nullptr;
	delete[] m_dataBuf;
	m_dataBuf = nullptr;
	delete[] m_init_S;
	m_init_S = nullptr;
	delete[] m_init_S_min;
	m_init_S_min = nullptr;
	delete[] m_init_S_tmp;
	m_init_S_tmp = nullptr;
	delete[] m_init_p1;
	m_init_p1 = nullptr;
	delete[] m_EN_cos;
	m_EN_cos = nullptr;
	delete[] m_EN_sin;
	m_EN_sin = nullptr;
	delete[] m_post_SNR;
	m_post_SNR = nullptr;
	delete[] m_pr_SNR;
	m_pr_SNR = nullptr;
	delete[] m_E_pr_SNR;
	m_E_pr_SNR = nullptr;
	delete[] m_v;
	m_v = nullptr;
	delete[] m_Gh1;
	m_Gh1 = nullptr;
	delete[] m_q;
	m_q = nullptr;
	delete[] m_pp;
	m_pp = nullptr;
	delete[] m_G;
	m_G = nullptr;
	delete[] m_M;
	m_M = nullptr;
	delete[] m_ns_storage;
	m_ns_storage = nullptr;
	delete[] m_cosen_local;
	m_cosen_local = nullptr;
	delete[] m_cosen_global;
	m_cosen_global = nullptr;
	delete[] m_plocal;
	m_plocal = nullptr;
	delete[] m_pglobal;
	m_pglobal = nullptr;
	delete[] m_h_global;
	m_h_global = nullptr;
	delete[] m_cosen;
	m_cosen = nullptr;
	delete[] m_old_cosen;
	m_old_cosen = nullptr;
	delete[] m_lamda_d;
	m_lamda_d = nullptr;
	delete[] m_integra;
	m_integra = nullptr;
}


LSA_denoise::LSA_denoise()
{
}

void LSA_denoise::Initialize(short wlen)
{
	m_lwlen = wlen;
	m_linc13 = m_lwlen/3;
	m_linc23 = m_linc13<<1;
	m_winData = new Complex_num[m_lwlen<<1];
	m_ns_hn = new int[m_lwlen];
	Gc.Initialize(m_linc23);
	MyN_fft.initial(m_linc23);

	for (int i = 1; i < (m_linc23+1); ++i) { 
		m_ns_hn[i - 1] = sqrt((0.5 - 0.5 * cos(2.0 * Pi*(i) / (m_linc23 + 1))) * 1073741824);
	}
}


void  LSA_denoise::Denoise_process( short* data_in, short* data_out , int blockInd)
{
	if (data_in == NULL)
		mcra_fatal("LSA_denoise::Denoise_process: null data_in");
	if (data_out == NULL)
		mcra_fatal("LSA_denoise::Denoise_process: null data_out");
	
	for (int i = 0; i < m_linc23; i++){
		m_winData[i].real = ((int64_t)data_in[i] *m_ns_hn[i]) >> 9;
		m_winData[i].imag = ((int64_t)data_in[i + m_linc13] * m_ns_hn[i]) >> 9;
	} 
 	MyN_fft.base4_fft(m_winData, 1);

	Gc.G_calculate_process(m_winData, blockInd);
	Gc.G_calculate_process(m_winData + m_linc23, blockInd+1);

	MyN_fft.base4_fft(m_winData, -1);
	MyN_fft.base4_fft(m_winData + m_linc23, -1);

	for (int i = 0; i < m_linc23; i++) {
		data_out[i]+= ((int64_t)m_winData[i].real * m_ns_hn[i]) >> 21;
		data_out[i+m_linc13] += ((int64_t)m_winData[i + m_linc23].real * m_ns_hn[i]) >> 21;
	}
}
LSA_denoise::~LSA_denoise()
{
	delete[] m_ns_hn;
	m_ns_hn = nullptr;
	delete[] m_winData;
	m_winData = nullptr;
}

Datablock_Read::Datablock_Read(int sample_rate,short channels,int MaxDataLen)
{
	m_channels = channels;
	m_sample_rate = sample_rate;
	Initial(MaxDataLen);
}

void Datablock_Read::Initial(int MaxDataLen) {
	m_maxdata = MaxDataLen;
	if (m_sample_rate < 23000)
	{
		m_wlen = 256;
		m_inc = 128;
	}
	else {
		m_wlen = 1024;
		m_inc = 512;
	}
	m_wlen15 = m_wlen +(m_wlen>>1);
	m_inc2 = m_inc <<1;

	m_inc_move = log(m_inc) / log(2);

	m_DoubDataBuffer = new short[MaxDataLen + m_wlen];
	m_data_in = new short[m_wlen15];
	m_data_storage = new short[m_wlen15];
	m_process_storage=new int[m_wlen15];
	m_buffer = new short[m_wlen15];
	m_data_out = new short[m_wlen15];
	m_data_resize = new short[MaxDataLen + m_wlen];

	memset(m_process_storage, 0, m_wlen15 * sizeof(int));
	memset(m_DoubDataBuffer, 0, (MaxDataLen + m_wlen) * sizeof(short));
	m_data_rest_length = 0;
	m_blockInd = 0;
	LSA.Initialize(m_wlen15);
}

void Datablock_Read::Data_procese(short* pInBuffer, short* pOutBuffer,int read_length,int& Out_Length){
	if (pInBuffer == NULL)
		mcra_fatal("Datablock_Read::Data_procese: null pInBuffer");
	if (pOutBuffer == NULL)
		mcra_fatal("Datablock_Read::Data_procese: null pOutBuffer");
	if (m_channels == 2)
		read_length = (read_length >> 1);
	int current_total_data = read_length + m_data_rest_length;

	if (current_total_data < m_wlen15) {
		if(m_channels==1)
		memcpy(m_data_storage + m_data_rest_length, pInBuffer, (read_length) * sizeof(short));
		else if (m_channels == 2) {
			for (int i = 0; i < read_length; i++) 
				m_data_storage[i + m_data_rest_length] = (pInBuffer[i << 1] + pInBuffer[(i << 1) + 1]) >> 1;
		}

		Out_Length = 0;
		m_data_rest_length = current_total_data;
	}
	else {
		memcpy(m_data_resize, m_data_storage, (m_data_rest_length) * sizeof(short));
		if (m_channels == 1){
			memcpy(m_data_resize + m_data_rest_length, pInBuffer, read_length * sizeof(short));
		}
		else if (m_channels == 2) {
			for (int i = 0; i < read_length; i++) {
				m_data_resize[i + m_data_rest_length] = (pInBuffer[i << 1] + pInBuffer[(i << 1) + 1]) >> 1;
			}
		}

		int data_use = 0;

		int current_frame = ((current_total_data - m_wlen) >> m_inc_move) + 1;

		current_frame = (current_frame | 1) - 1;

		m_data_rest_length = current_total_data -( current_frame*m_inc);
		Out_Length = current_frame * m_inc;


		for (int k = 0; k < current_frame-1; k+=2, m_blockInd+=2){
			memcpy(m_data_in, m_data_resize + data_use, m_wlen15 * sizeof(short));
			memset(m_data_out, 0, sizeof(short)*m_wlen15);
			LSA.Denoise_process(m_data_in, m_data_out, m_blockInd);
			int block_inc = k * m_inc;
			if (current_frame == 2) {
				for (int i = 0; i < m_inc; i++) {
					m_DoubDataBuffer[block_inc + i] += m_data_out[i];
					m_DoubDataBuffer[block_inc + i] += m_process_storage[i];
					m_DoubDataBuffer[block_inc + m_inc+i] += m_data_out[i+m_inc];
					m_process_storage[i] = m_data_out[i+m_wlen];
				}
			}
			else {
				if (k == 0) {                    
					for (int i = 0; i < m_inc; i++) {
						m_DoubDataBuffer[block_inc + i] += m_data_out[i];
						m_DoubDataBuffer[block_inc + i] += m_process_storage[i];
						m_DoubDataBuffer[block_inc + m_inc + i] += m_data_out[i + m_inc];
						m_DoubDataBuffer[block_inc + m_wlen + i ] += m_data_out[i + m_wlen];
					}
				}
				else if (k == current_frame-2 ){
					for (int i = 0; i < m_inc; i++) {   
						m_DoubDataBuffer[block_inc + i] += m_data_out[i];
						m_DoubDataBuffer[block_inc + m_inc + i] += m_data_out[i + m_inc];
						m_process_storage[i] = m_data_out[i + m_wlen];
					}
				}
				else {
					for (int i = 0; i < m_wlen15; i++) {
						m_DoubDataBuffer[block_inc + i] += m_data_out[i];
					}
				}
			}
			if (m_channels == 1) {
				for (int i = 0; i < Out_Length; i++) {
					pOutBuffer[i] = m_DoubDataBuffer[i];
				}
			}
			else if (m_channels == 2){
				for (int i = 0; i < Out_Length; i++) {
					pOutBuffer[i << 1] = m_DoubDataBuffer[i];
					pOutBuffer[(i << 1) + 1] = m_DoubDataBuffer[i];
				}
			}
			data_use += m_wlen;
		}
		if (m_channels == 2)
			Out_Length = Out_Length << 1;  
		memcpy(m_data_storage, m_data_resize + data_use, (m_data_rest_length) * sizeof(short));
		memset(m_DoubDataBuffer, 0, (m_maxdata + m_wlen) * sizeof(short));
	}
}

	 
Datablock_Read::~Datablock_Read()
{
	delete[] m_data_storage;
	m_data_storage = nullptr;
	delete[] m_data_in;
	m_data_in = nullptr;
	delete[] m_data_resize;
	m_data_resize = nullptr;
	delete[] m_buffer;
	m_buffer = nullptr;
	delete[] m_process_storage;
	m_process_storage = nullptr;
	delete[] m_DoubDataBuffer;
	m_DoubDataBuffer = nullptr;
	delete[] m_data_out;
	m_data_out = nullptr;
}



void procese(const char *szFileIn, const char *szFileOut,int sample_rate)
{
	FILE *fpIn=fopen(szFileIn, "rb");
	if (NULL == fpIn)
		mcra_fatal("procese: cannot open input file");

	FILE *fpOut=fopen(szFileOut, "wb"); 
	if (NULL == fpOut)
		mcra_fatal("procese: cannot open output file");

	short pInBuffer[section_max];
	memset(pInBuffer, 0, section_max*sizeof(short));
	short pOutBuffer[section_max+frame_max];
	memset(pOutBuffer, 0, (section_max + frame_max) * sizeof(short));

	int Out_Length = 0;	
	short voice_channel = 1;
	int MaxDataLen = 10000;
	Datablock_Read dtr(sample_rate, voice_channel, MaxDataLen);

	clock_t start = clock();
	while (!feof(fpIn))
	{
		int read_length = 160; 
		int nread = fread(pInBuffer, sizeof(short), read_length, fpIn);
		if (nread > 0)
		{
			dtr.Data_procese(pInBuffer, pOutBuffer, read_length, Out_Length);
			if (Out_Length > 1) {
				fwrite(pOutBuffer, sizeof(short), Out_Length, fpOut);
			}
			memset(pInBuffer, 0, MaxDataLen * sizeof(short));
			memset(pOutBuffer, 0, (MaxDataLen + frame_max) * sizeof(short));
		}
	}
	clock_t finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("%f seconds", duration);
	fclose(fpIn);
	fclose(fpOut);
}

int main(const int argc, const char * argv[]) {
	procese("60db_16k.pcm", "60db_16k_out.pcm",16000);
	for (int i = 0; i < argc; i++) {
		printf("%s\n", argv[i]);
	}
	printf("降噪完成\n");
	return 0;
}